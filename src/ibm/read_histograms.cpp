// file to read in data
// and to produce histograms of the data

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <vector>
#include "auxiliary.h"
#include <gsl/gsl_histogram.h>

using namespace std;

const size_t number_columns = 2;

template <class Container>

// function to split filenames into folders
void split_folders(const std::string& str, Container &cont, char delim='/')
{
    stringstream ss(str);
    string token;

    // ok, split the string according to the delimiter
    while (getline(ss, token, delim))
    {
            cont.push_back(token);
    }
}


// opens the file with all the values
void initFile(
        int argc, 
        char **argv, 
        ifstream &file, 
        string &base_path)
{
    // first get the path name of the file
    string filename(argv[1]);

    // then split the path name into all folders, apart from the last
    vector <string> folders;
    split_folders(filename, folders, '/');

    // now iterate over folders to get everything but the last 
    base_path = "";

    for (size_t i = 0; i < folders.size() - 1; ++i)
    {
        base_path += folders[i] + "/";
    }

    file.open(argv[1]);
}

// finds the minima and the maxima per specific column
void findExtremes(double *max, 
        double *min, 
        ifstream &file, 
        unsigned long &file_begin,
        string &header)
{
    // go through the number of columns and 
    // calculate minima and maxima
    for (size_t i = 0; i < number_columns; ++i)
    {
        max[i] = 0;
        min[i] = 0;
    }

    string line;

    // store the beginning of the file
    file_begin = file.tellg();

    // get the first line containing the headers
    std::getline(file, line);

    // store the file header
    header = line;

    // read the lines
    while (std::getline(file, line))
    {
        // read each line into memory
        // since we need to split it
        // into csv values
        istringstream linestream(line);

        size_t itemnum = 0;

        string item;

        // do the splitting
        while (std::getline(linestream, item, ';'))
        {
            if (itemnum > number_columns)
            {
                break;
            }
   
            // skip first column 
            // (which contains the generation number)
            if (itemnum > 0)
            {
                // obtain the value
                double curval = atof(item.c_str());

                // check whether this value is within range
                assert((itemnum - 1) < number_columns);

                // update maxima and minima
                curval > max[itemnum - 1] ? 
                    max[itemnum - 1] = curval 
                    : 
                    (min[itemnum - 1] > curval ? 
                        min[itemnum - 1] = curval : 0);
            }

            ++itemnum;
        }
    }

    file.clear();
}

void initHistograms(gsl_histogram **histograms, 
        size_t histogram_max, 
        double *max, 
        double *min)
{
    for (size_t i = 0; i < histogram_max; ++i)
    {
        histograms[i] = gsl_histogram_alloc(500);
        cout << "max " << i << ": " << max[i] << endl;
        cout << "min" << i << ": " << min[i] << endl;

        if (min[i] == max[i])
        {
            max[i] += 1;
        }

        // add 10% to the outer boundarys
        min[i] -= 0.1 * fabs(min[i]);
        max[i] += 0.1 * fabs(max[i]);

        // set the ranges of the histogram
        gsl_histogram_set_ranges_uniform(
                histograms[i], 
                min[i], 
                max[i]);
    }
}

// reset the histogram
void resetHistograms(
        gsl_histogram **histograms, 
        size_t histogram_max)
{
    for (size_t i = 0; i < histogram_max; ++i)
    {
        gsl_histogram_reset(histograms[i]);
    }
}

// write histograms to output file
void writeHistograms(
        gsl_histogram **histograms, 
        size_t histogram_max, 
        size_t generation, 
        FILE *file,
        string file_header)
{
    stringstream header_line(file_header);

    string header_item;

    size_t itemnum = 0;

    stringstream histostrings[number_columns];

    // loop through the lines 
    while (std::getline(header_line, header_item, ';')) 
    {
        if (itemnum > number_columns)
        {
            break;
        }

        if (itemnum > 0)
        {
            histostrings[itemnum-1] << itos(generation) << ";" 
                << header_item << ";%f;";
        }

        ++itemnum;
    }

    for (size_t i = 0; i < histogram_max; ++i)
    {
        gsl_histogram_fprintf(
                file, 
                histograms[i],"%e;",histostrings[i].str().c_str());
    }
}

// fill the histograms
void fillHistograms(
        double * max, 
        double * min, 
        ifstream &file, 
        unsigned long file_begin,
        string &file_header,
        string &hist_file_name
        )
{
    gsl_histogram * histograms[number_columns];

    FILE *myfile;
    
    if ( ! (myfile = fopen(hist_file_name.c_str(),"w")))
    {
        cout << "cannot open histograms for writing!" << endl;
        exit(1);
    }

    // initialize the histograms
    initHistograms(histograms, number_columns, max, min);

    string line;

    int generation = 0;

    file.seekg(file_begin, ios::beg);

    // skip the first line
    std::getline(file, line);

    // loop through lines
    while (std::getline(file, line))
    {
        // read each line into memory
        // since we need to split it
        // into csv values
        istringstream linestream(line);

        size_t itemnum = 0;

        string item;
       
        // split elements of each line
        while (std::getline(linestream, item, ';'))
        {
            // check whether the generation is still the same
            // if not, write the histograms of the previous 
            // generation to the file 
            // and reset histograms
            if (itemnum == 0 && generation != atoi(item.c_str()))
            {
                cout << "write generation: " << generation << endl;

                // write the histograms to the output file
                writeHistograms(histograms, 
                        number_columns, 
                        generation,
                        myfile,
                        file_header
                        );
                
                resetHistograms(histograms, number_columns);
                
                generation = atoi(item.c_str());
            }
            else if (itemnum > 0)
            {
                // increment the corresponding histogram
                gsl_histogram_increment(
                        histograms[itemnum - 1], 
                        atof(item.c_str())
                        );
            }
            ++itemnum;
        }
    }
    
    writeHistograms(histograms, 
            number_columns,
            generation,
            myfile,
            file_header);
}

// the guts of the code
int main(int argc, char **argv)
{
    ifstream file;

    string base_path;

    // open the file
    initFile(argc, argv, file, base_path);
    
    base_path += "histograms.csv";

    // minima and maxima of each column
    double max[number_columns];
    double min[number_columns];

    // store the starting point of the file
    unsigned long file_begin;

    // store the file header which can be later 
    // written to the histogram output file
    string file_header;

    findExtremes(max, min, file, file_begin, file_header);

    cout << file_begin << endl;

    fillHistograms(max, min, file, file_begin, file_header, base_path);
}
