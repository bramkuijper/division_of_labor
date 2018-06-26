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
// also stores the folder in which the file is present
// as we need to write the output file to that folder as well
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
        string &header,
        size_t const number_columns
        )
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

    // process the header line
    // if the first character is a digit
    // then this not a header line
    if (isdigit(line[0]))
    {
        // hence we have to make a header ourselves
        // first item is generation
        header = "generation";

        // next items are the individual traits, which
        // we just label trait0, trait1, trait2, ...
        for (size_t col_i = 0; col_i < number_columns; ++col_i)
        {
            stringstream ss;

            ss << ";trait" << col_i;
            header += ss.str();
        }
    }
    else
    {
        // yes, header is present, so use it
        header = line;
    }


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
        size_t number_columns, 
        double *max, 
        double *min)
{
    for (size_t i = 0; i < number_columns; ++i)
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
        size_t number_columns, 
        size_t generation, 
        FILE *file,
        string file_header
        )
{
    // the output of the histogram is
    // bin[i-1];bin[i];generation;variable_name;count;

    // make a stringstream to do some splitting work
    // on the header line
    stringstream header_line(file_header);

    string header_item;

    size_t itemnum = 0;

    // store a bunch of string streams 
    // that each represents how output should 
    // work for each histogram
    stringstream histostrings[number_columns];

    // loop through the header items
    // item by item (split by ";")
    while (std::getline(header_line, header_item, ';')) 
    {
        if (itemnum > number_columns)
        {
            break;
        }

        // now store the trait name 
        if (itemnum > 0) // skip the first item (as this is generation)
        {
            // now generate the bin format for the histogram
            // this should be 
            histostrings[itemnum-1] << itos(generation) << ";" 
                << header_item << ";%f";
        }

        ++itemnum;
    }

    cout << histostrings[2].str() << endl;

    // now print the histogram
    for (size_t i = 0; i < number_columns; ++i)
    {
        gsl_histogram_fprintf(
                file, // stream
                histograms[i], // the histogram object
                "%f;", //the range format
                histostrings[i].str().c_str() // the bin format);
                ); 
    }
}

// fill the histograms
void fillHistograms(
        double * max, 
        double * min, 
        ifstream &file, 
        unsigned long file_begin,
        string &file_header,
        string &hist_file_name,
        size_t const number_columns
        )
{
    gsl_histogram * histograms[number_columns];


    // file to write output to
    FILE *myfile;
    
    if ( ! (myfile = fopen(hist_file_name.c_str(),"w")))
    {
        cout << "cannot open histograms for writing!" << endl;
        exit(1);
    }

    // make header for the histogram 
    string output_file_header = "bin_start;bin_end;generation;traitname;count\n";

    // write the header to the histogram
    fprintf(myfile, "%s", output_file_header.c_str());

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
            else if (itemnum > 0 && itemnum < number_columns + 1)
            {
                cout << "\"" << atof(item.c_str()) << "\"" << endl;
                // increment the corresponding histogram
                gsl_histogram_increment(
                        histograms[itemnum - 1], 
                        atof(item.c_str())
                        );

                cout << "done" << endl;
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

    // variable to store the folder name 
    // in which we find the histogram 
    // file. We need to store the output file there
    string base_path;

    // get number of columns present in the file
    // from the command line
    size_t number_columns = atoi(argv[2]);

    string output_file = argv[3];

    // open the file
    // and also stores the basename of the file
    initFile(argc, argv, file, base_path);
    
    // add the name of the output file
    base_path = output_file;

    // minima and maxima of each column
    double max[number_columns];
    double min[number_columns];

    // store the starting point of the file
    unsigned long file_begin;

    // store the file header which can be later 
    // written to the histogram output file
    string file_header;

    findExtremes(max, 
            min, 
            file, 
            file_begin, 
            file_header, 
            number_columns);

    fillHistograms(max, 
            min, 
            file, 
            file_begin, 
            file_header, 
            base_path,
            number_columns);

    file.close();

    return(0);
}
