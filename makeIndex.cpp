#include <ctime>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include "GridIndex.h"

#define LATITUDE_MIN -90.0
#define LATITUDE_MAX 90.0
#define LONGITUDE_MIN -176.3
#define LONGITUDE_MAX 177.5

/*
    ./makeIndex Gowalla_totalCheckins.txt index_10.txt pp_data.txt 10
    ./makeIndex Gowalla_totalCheckins.txt index_50.txt pp_data.txt 50
    ./makeIndex Gowalla_totalCheckins.txt index_100.txt pp_data.txt 100
    ./makeIndex Gowalla_totalCheckins.txt index_150.txt pp_data.txt 150
    ./makeIndex Gowalla_totalCheckins.txt index_200.txt pp_data.txt 200
*/

using namespace std;

int duplicate_elimination(char* data_path, char* data_path_new) {
    // read the original dataset from data_path
    // eliminate duplicates by deleting the corresponding lines
    // write the dataset without duplicates into data_path_new
    fstream data_file;
    data_file.open(data_path, ios::in);
    if (!data_file.is_open()) {
        cout << "FAILED TO OPEN " << data_path << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    cout << data_path << " opened successfully." << endl;
    string raw_spatial_point;
    list<SpatialPoint> spatial_data;
    cout << "Loading raw data from " << data_path << "." << endl;
    while (getline(data_file, raw_spatial_point)) {
        istringstream raw_sp_sstream(raw_spatial_point);
        string temp;
        double latitude, longitude;
        int location_id;
        // ignore the first two columns
        raw_sp_sstream >> temp;
        raw_sp_sstream >> temp;
        raw_sp_sstream >> latitude >> longitude >> location_id;
        if (latitude < LATITUDE_MIN || latitude > LATITUDE_MAX ||
            longitude < LONGITUDE_MIN || longitude > LONGITUDE_MAX)
            continue;
        spatial_data.push_back(SpatialPoint(latitude, longitude, location_id));
    }
    data_file.close();
    // sorting spatial points in the ascending order of latitudes, then longitudes, and finally location_ids
    spatial_data.sort(compare_spatial_point);
    fstream data_file_new;
    data_file_new.open(data_path_new, ios::out);
    if (!data_file_new.is_open()) {
        cout << "FAILED TO OPEN " << data_path_new << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    cout << "Removing duplications and writing preprocessed data to " << data_path_new << "." << endl;
    data_file_new.precision(PRECISION);
    data_file_new << fixed;
    auto it_i = spatial_data.begin();
    while (!spatial_data.empty()) {
        auto it_j = it_i;
        while (next(it_j)->get_latitude() == it_j->get_latitude() &&
               next(it_j)->get_longitude() == it_j->get_longitude() &&
               next(it_j) != spatial_data.end())
            it_j++;
        data_file_new << it_i->get_latitude() << "\t"
                      << it_i->get_longitude() << "\t"
                      << it_i->get_location_id() << endl;
        spatial_data.erase(it_i, next(it_j));
        it_i = spatial_data.begin();
    }
    data_file_new.close();
    cout << "Data preprocessing completed successfully." << endl;
    return 0;
}

string create_index(char* data_path_new, char* index_path, int n) {
    // To create a grid index and save it to file on "index_path".
    // The output file should contain exactly n*n lines. If there is no point in the cell, just leave it empty after ":".
    fstream data_file_new;
    data_file_new.open(data_path_new, ios::in);
    if (!data_file_new.is_open()) {
        cout << "FAILED TO OPEN " << data_path_new << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    cout << data_path_new << " opened successfully." << endl;
    string raw_spatial_point;
    list<SpatialPoint> spatial_data;
    cout << "Loading preprocessed data from " << data_path_new << "." << endl;
    while (getline(data_file_new, raw_spatial_point)) {
        istringstream raw_sp_sstream(raw_spatial_point);
        double latitude, longitude;
        int location_id;
        raw_sp_sstream >> latitude >> longitude >> location_id;
        spatial_data.push_back(SpatialPoint(latitude, longitude, location_id));
    }
    data_file_new.close();
    fstream index_file;
    index_file.open(index_path, ios::out);
    if (!index_file.is_open()) {
        cout << "FAILED TO OPEN " << index_path << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    index_file.precision(PRECISION);
    index_file << fixed;
    double cell_range_x = (LATITUDE_MAX - LATITUDE_MIN) / n;
    double cell_range_y = (LONGITUDE_MAX - LONGITUDE_MIN) / n;
    index_file << LATITUDE_MIN << " " << LATITUDE_MAX << " " << LONGITUDE_MIN << " " << LONGITUDE_MAX << endl;
    cout << "Creating the index and writing it to " << index_path << "." << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double x_min = LATITUDE_MIN + cell_range_x * j;
            double x_max = LATITUDE_MIN + cell_range_x * (j + 1);
            double y_min = LONGITUDE_MIN + cell_range_y * i;
            double y_max = LONGITUDE_MIN + cell_range_y * (i + 1);
            index_file << "Cell " << i << "," << j << ":";
            bool top_border_cell = false;
            bool left_border_cell = false;
            if (i == 0)
                top_border_cell = true;
            if (j == 0)
                left_border_cell = true;
            auto it = spatial_data.begin();
            bool last_point_considered = false;
            while(!spatial_data.empty() && !last_point_considered) {
                bool point_in_cell = false;
                if (left_border_cell) {
                    if ((it->get_latitude() > x_min || is_equal(it->get_latitude(), x_min)) &&
                        (it->get_latitude() < x_max || is_equal(it->get_latitude(), x_max))) {
                        if (top_border_cell) {
                            if ((it->get_longitude() > y_min || is_equal(it->get_longitude(), y_min)) &&
                                (it->get_longitude() < y_max || is_equal(it->get_longitude(), y_max)))
                                point_in_cell = true;
                        }
                        else {
                            if (it->get_longitude() > y_min &&
                                (it->get_longitude() < y_max || is_equal(it->get_longitude(), y_max)))
                                point_in_cell = true;
                        }
                    }
                }
                else {
                    if (it->get_latitude() > x_min &&
                        (it->get_latitude() < x_max || is_equal(it->get_latitude(), x_max))) {
                        if (top_border_cell) {
                            if ((it->get_longitude() > y_min || is_equal(it->get_longitude(), y_min)) &&
                                (it->get_longitude() < y_max || is_equal(it->get_longitude(), y_max)))
                                point_in_cell = true;
                        }
                        else {
                            if (it->get_longitude() > y_min &&
                                (it->get_longitude() < y_max || is_equal(it->get_longitude(), y_max)))
                                point_in_cell = true;
                        }
                    }
                }
                if (it == prev(spatial_data.end()))
                    last_point_considered = true;
                else if (next(it)->get_latitude() > x_max)
                    last_point_considered = true;
                if (point_in_cell) {
                    index_file << " " << it->get_location_id()
                               << "_" << it->get_latitude()
                               << "_" << it->get_longitude();
                    it = spatial_data.erase(it);
                    it--;
                }
                it++;
            }
            index_file << endl;
        }
    }
    index_file.close();
    return "";
}

// this main() is for quickly creating index files for n = 10, 50, 100, 150, 200
// (assume data_path, index_path, and data_path_new are known.)

/*
int main() {
    char data_path[] = "Gowalla_totalCheckins.txt";
    char data_path_new[] = "pp_data.txt";
    duplicate_elimination(data_path, data_path_new);
    vector<int> index_n {10, 50, 100, 150, 200};
    for (auto it = index_n.begin(); it != index_n.end(); it++) {
        clock_t s, t;
        s = clock();
        char index_path[20];
        sprintf(index_path, "index_%d.txt", *it);
        create_index(data_path_new, index_path, *it);
        t = clock();
        cout << fixed << "index construction time: " << (double)(t-s) << endl;
    }
    return 0;
}
*/

int main(int argc, char** argv) {
    if (argc != 5){
        cout << "Ussage: " << argv[0] << " DATA_PATH INDEX_PATH DATA_PATH_NEW N" << endl;
        
        // DATA_PATH(char *): the file path of Gowalla_totalCheckins.txt
        // INDEX_PATH(char *): the output file path of the grid index
        // DATA_PATH_NEW(char *): the file path of the dataset without duplicates
        // N(integer): the grid index size
        
        return -1;
    }
    duplicate_elimination(argv[1], argv[3]);
    clock_t s, t;
    s = clock();
    create_index(argv[3], argv[2], atoi(argv[4]));
    t = clock();
    cout << fixed << "index construction time: " << (double)(t-s) << endl;
    return 0;
}
