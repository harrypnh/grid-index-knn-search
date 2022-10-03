#include <ctime>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <cfloat>
#include <random>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "GridIndex.h"

// change this const to experiment different M values for the R-tree
#define RTREE_MAX_ELEMENTS 4096
#define BUFFER_RATE 0.7

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::d2::point_xy<double> Point;
typedef pair<Point, int> RTreePoint;
typedef bgi::rtree<RTreePoint, bgi::rstar<RTREE_MAX_ELEMENTS>> RTree;

string knn_grid(double x, double y, char* index_path, int k, int n) {
    // to get the 5-NN result with the help of the grid index
    // Please store the 5-NN results by a String of location ids, like "11, 789, 125, 2, 771"
    GridIndex grid_index(index_path, n);
    return grid_index.knn_search(x, y, k);
}

string knn_grid_fast(double x, double y, char* index_path, int k, int n) {
    // to get the 5-NN result with the help of the grid index
    // Please store the 5-NN results by a String of location ids, like "11, 789, 125, 2, 771"
    GridIndex grid_index(index_path, n);
    return grid_index.knn_search_fast(x, y, k);
}

string knn_grid_compress(double x, double y, char* index_path, int k, int n) {
    // to get the 5-NN result with the help of the grid index
    // Please store the 5-NN results by a String of location ids, like "11, 789, 125, 2, 771"
    GridIndexCompress grid_index_compress(index_path, n);
    return grid_index_compress.knn_search(x, y, k);
}

string knn_grid_disk(double x, double y, char* index_path, int k, int n, double buffer_rate) {
    // to get the 5-NN result with the help of the grid index
    // Please store the 5-NN results by a String of location ids, like "11, 789, 125, 2, 771"
    GridIndexDisk grid_index_disk(index_path, n, buffer_rate);
    return grid_index_disk.knn_search(x, y, k);
}

RTree preload_rtree(char* data_path_new) {
    fstream data_file_new;
    data_file_new.open(data_path_new, ios::in);
    if (!data_file_new.is_open()) {
        cout << "FAILED TO OPEN " << data_path_new << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    string raw_spatial_point;
    vector<RTreePoint> spatial_data;
    while (getline(data_file_new, raw_spatial_point)) {
        istringstream raw_sp_sstream(raw_spatial_point);
        double latitude, longitude;
        int location_id;
        raw_sp_sstream >> latitude >> longitude >> location_id;
        Point spatial_point(latitude, longitude);
        spatial_data.push_back(make_pair(spatial_point, location_id));
    }
    data_file_new.close();
    RTree rtree(spatial_data.begin(), spatial_data.end());
    return rtree;
}

string knn_rtree_preload(RTree rtree, double x, double y, int k) {
    list<RTreePoint> knn;
    Point query_point(x, y);
    for (auto it = rtree.qbegin(bgi::nearest(query_point, k)); it != rtree.qend(); it++)
        knn.push_back(*it);
    ostringstream results_ss;
    for (auto it = knn.begin(); it != knn.end(); it++) {
        if (it != prev(knn.end()))
            results_ss << it->second << ", ";
        else
            results_ss << it->second;
    }
    return results_ss.str();
}

string knn_rtree(double x, double y, char* data_path_new, int k) {
    RTree rtree = preload_rtree(data_path_new);
    return knn_rtree_preload(rtree, x, y, k);
}

list<SpatialPoint> preload_dataset(char* data_path_new) {
    fstream data_file_new;
    data_file_new.open(data_path_new, ios::in);
    if (!data_file_new.is_open()) {
        cout << "FAILED TO OPEN " << data_path_new << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    string raw_spatial_point;
    list<SpatialPoint> spatial_data;
    while (getline(data_file_new, raw_spatial_point)) {
        istringstream raw_sp_sstream(raw_spatial_point);
        double latitude, longitude;
        int location_id;
        raw_sp_sstream >> latitude >> longitude >> location_id;
        spatial_data.push_back(SpatialPoint(latitude, longitude, location_id));
    }
    data_file_new.close();
    return spatial_data;
}

string knn_linear_scan_preload(list<SpatialPoint> spatial_data, double x, double y, int k) {
    SpatialPoint query_point(x, y, 0, 0);
    list<SpatialPoint> knn;
    double knn_max_distance = DBL_MAX;
    for (auto point_it = spatial_data.begin(); point_it != spatial_data.end(); point_it++) {
        point_it->calculate_query_distance(query_point);
        if (knn.size() < k) {
            knn.push_back(SpatialPoint(point_it->get_latitude(), point_it->get_longitude(), point_it->get_location_id(), point_it->get_query_distance()));
            if (knn.size() == 1 || (knn.size() > 1 && knn_max_distance < point_it->get_query_distance()))
                knn_max_distance = point_it->get_query_distance();
        }
        else if (point_it->get_query_distance() < knn_max_distance) {
            auto farthest_nn_it = knn.begin();
            for (auto knn_it = knn.begin(); knn_it != knn.end(); knn_it++)
                if (knn_it->get_query_distance() > farthest_nn_it->get_query_distance())
                    farthest_nn_it = knn_it;
            knn.erase(farthest_nn_it);
            knn.push_back(SpatialPoint(point_it->get_latitude(), point_it->get_longitude(), point_it->get_location_id(), point_it->get_query_distance()));
            farthest_nn_it = knn.begin();
            for (auto knn_it = knn.begin(); knn_it != knn.end(); knn_it++)
                if (knn_it->get_query_distance() > farthest_nn_it->get_query_distance())
                    farthest_nn_it = knn_it;
            knn_max_distance = farthest_nn_it->get_query_distance();
        }
    }
    knn.sort(compare_query_distance);
    ostringstream results_ss;
    for (auto it = knn.begin(); it != knn.end(); it++) {
        if (it != prev(knn.end()))
            results_ss << it->get_location_id() << ", ";
        else
            results_ss << it->get_location_id();
    }
    return results_ss.str();
}

string knn_linear_scan(double x, double y, char* data_path_new, int k) {
    // to get the 5-NN result by linear scan
    // Please store the 5-NN results by a String of location ids, like "11, 789, 125, 2, 771"
    list<SpatialPoint> spatial_data = preload_dataset(data_path_new);
    return knn_linear_scan_preload(spatial_data, x, y, k);
}

// this main() is for evaluate the running time of all knn search methods on random query points
// (assume data_path_new and all index_path are known.)
/*
int main(int argc, char** argv) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " RANDOM_QUERY_SIZE EVAL_DATA_PATH K_MAX" << endl;
    
        // RANDOM_QUERY_SIZE(integer): the number of random query points to generate
        // EVAL_DATA_PATH(char *): the file path to output statistics data for evaluation
        // K_MAX(integer): the k_max value for k-NN search for k from 1 to k_max
    
        return -1;
    }
    const int RANDOM_QUERY_SIZE = atoi(argv[1]);
    const int K_MAX = atoi(argv[3]);
    char* eval_data_path = argv[2];
    cout << fixed;
    char data_path_new[] = "pp_data.txt";
    vector<SpatialPoint> random_query;
    random_query.reserve(RANDOM_QUERY_SIZE);
    uniform_real_distribution<double> x_range(-90.0, 90.0);
    uniform_real_distribution<double> y_range(-176.3, 177.5);
    mt19937 random_generator;
    random_generator.seed(random_device{}());
    for (int i = 0; i < RANDOM_QUERY_SIZE; i++) {
        double x_q = x_range(random_generator);
        double y_q = y_range(random_generator);
        SpatialPoint q(x_q, y_q, 0, 0);
        random_query.push_back(q);
    }
    vector<int> index_n {10, 50, 100, 150, 200};
    vector<GridIndex> grid_indexes;
    vector<GridIndexCompress> grid_indexes_compress;
    // vector<GridIndexDisk> grid_indexes_disk;
    grid_indexes.reserve(index_n.size());
    grid_indexes_compress.reserve(index_n.size());
    // grid_indexes_disk.reserve(index_n.size());
    
    cout << "Loading the preprocessed data" << endl;
    list<SpatialPoint> spatial_data = preload_dataset(data_path_new);
    cout << endl;

    cout << "Loading the grid index files (normal mode)" << endl;
    vector<double> grid_index_loading_time;
    grid_index_loading_time.reserve(index_n.size());
    for (int i = 0; i < index_n.size(); i++) {
        char index_path[20];
        sprintf(index_path, "index_%d.txt", index_n[i]);
        cout << "    Loading " << index_path << endl;
        clock_t s = clock();
        GridIndex grid_index(index_path, index_n[i]);
        clock_t t = clock();
        grid_index_loading_time.push_back(double(t - s));
        grid_indexes.push_back(grid_index);
    }
    cout << endl;
    cout << "Loading the grid index files (compress mode)" << endl;
    vector<double> grid_index_compress_loading_time;
    grid_index_compress_loading_time.reserve(index_n.size());
    for (int i = 0; i < index_n.size(); i++) {
        char index_path[20];
        sprintf(index_path, "index_%d.txt", index_n[i]);
        cout << "    Loading " << index_path << endl;
        clock_t s = clock();
        GridIndexCompress grid_index_compress(index_path, index_n[i]);
        clock_t t = clock();
        grid_indexes_compress.push_back(grid_index_compress);
        grid_index_compress_loading_time.push_back(double(t - s));
    }
    cout << endl;
    // cout << "Loading the grid index files (disk mode)" << endl;
    // vector<double> grid_index_disk_loading_time;
    // grid_index_disk_loading_time.reserve(index_n.size());
    // for (int i = 0; i < index_n.size(); i++) {
    //     char index_path[20];
    //     sprintf(index_path, "index_%d.txt", index_n[i]);
    //     cout << "    Loading " << index_path << endl;
    //     clock_t s = clock();
    //     GridIndexDisk grid_index_disk(index_path, index_n[i], BUFFER_RATE);
    //     clock_t t = clock();
    //     grid_indexes_disk.push_back(grid_index_disk);
    //     grid_index_disk_loading_time.push_back(double(t - s));
    // }
    // cout << endl;
    cout << "Loading the preprocessed dataset and building an r-tree" << endl;
    clock_t s = clock();
    RTree rtree = preload_rtree(data_path_new);
    clock_t t = clock();
    double rtree_loading_time = double(t - s);
    
    cout << endl;
    vector<double> knn_linear_scan_time;
    knn_linear_scan_time.reserve(K_MAX);
    vector<double> empty_time;
    empty_time.reserve(K_MAX);
    vector<vector<double>> knn_grid_time(index_n.size(), empty_time);
    vector<vector<double>> knn_grid_fast_time(index_n.size(), empty_time);
    vector<vector<double>> knn_grid_compress_time(index_n.size(), empty_time);
    // vector<vector<double>> knn_grid_disk_time(index_n.size(), empty_time);
    vector<double> knn_rtree_time;
    knn_rtree_time.reserve(K_MAX);

    vector<string> knn_linear_scan_results;
    knn_linear_scan_results.reserve(K_MAX * RANDOM_QUERY_SIZE);
    vector<string> empty_results;
    empty_results.reserve(K_MAX * RANDOM_QUERY_SIZE);
    vector<vector<string>> knn_grid_results(index_n.size(), empty_results);
    vector<vector<string>> knn_grid_fast_results(index_n.size(), empty_results);
    vector<vector<string>> knn_grid_compress_results(index_n.size(), empty_results);
    // vector<vector<string>> knn_grid_disk_results(index_n.size(), empty_results);
    vector<string> knn_rtree_results;
    knn_rtree_results.reserve(K_MAX * RANDOM_QUERY_SIZE);

    for (int k = 1; k <= K_MAX; k++) {
        cout << k << "-NN search on " << RANDOM_QUERY_SIZE << " random query points" << endl;
        double knn_linear_scan_time_avg = 0;
        cout << "    " << "using linear scan" << endl;
        for (int i = 0; i < random_query.size(); i++) {
            clock_t s = clock();
            string linear_result = knn_linear_scan_preload(spatial_data, random_query[i].get_latitude(), random_query[i].get_longitude(), k);
            clock_t t = clock();
            knn_linear_scan_time_avg += double(t - s) / random_query.size();
            knn_linear_scan_results.push_back(linear_result);
        }
        knn_linear_scan_time.push_back(knn_linear_scan_time_avg);
        
        cout << "    " << "using grid indexes (normal, fast, compress, disk) scan" << endl;
        for (int n_i = 0; n_i < index_n.size(); n_i++) {
            double knn_grid_time_n_avg = 0;
            double knn_grid_fast_time_n_avg = 0;
            double knn_grid_compress_time_n_avg = 0;
            // double knn_grid_disk_time_n_avg = 0;
            for (int i = 0; i < random_query.size(); i++) {
                char index_path[20];
                clock_t s, t;
                s = clock();
                string grid_result = grid_indexes[n_i].knn_search(random_query[i].get_latitude(), random_query[i].get_longitude(), k);
                t = clock(); 
                knn_grid_time_n_avg += double(t - s) / random_query.size();
                knn_grid_results[n_i].push_back(grid_result);
                s = clock();
                string grid_fast_result = grid_indexes[n_i].knn_search_fast(random_query[i].get_latitude(), random_query[i].get_longitude(), k);
                t = clock();
                knn_grid_fast_time_n_avg += double(t - s) / random_query.size();
                knn_grid_fast_results[n_i].push_back(grid_fast_result);
                s = clock();
                string grid_compress_result = grid_indexes_compress[n_i].knn_search(random_query[i].get_latitude(), random_query[i].get_longitude(), k);
                t = clock();
                knn_grid_compress_time_n_avg += double(t - s) / random_query.size();
                knn_grid_compress_results[n_i].push_back(grid_compress_result);
                // cout << random_query[i].get_latitude() << ", " << random_query[i].get_longitude() << endl;
                // s = clock();
                // string grid_disk_result = grid_indexes_disk[n_i].knn_search(random_query[i].get_latitude(), random_query[i].get_longitude(), k);
                // t = clock();
                // knn_grid_disk_time_n_avg += double(t - s) / random_query.size();
                // knn_grid_disk_results[n_i].push_back(grid_disk_result);
            }
            knn_grid_time[n_i].push_back(knn_grid_time_n_avg);
            knn_grid_fast_time[n_i].push_back(knn_grid_fast_time_n_avg);
            knn_grid_compress_time[n_i].push_back(knn_grid_compress_time_n_avg);
            // knn_grid_disk_time[n_i].push_back(knn_grid_disk_time_n_avg);
        }

        double knn_rtree_time_avg = 0;
        cout << "    " << "using rtree scan" << endl;
        for (int i = 0; i < random_query.size(); i++) {
            clock_t s = clock();
            string rtree_result = knn_rtree_preload(rtree, random_query[i].get_latitude(), random_query[i].get_longitude(), k);
            clock_t t = clock();
            knn_rtree_time_avg += double(t - s) / random_query.size();
            knn_rtree_results.push_back(rtree_result);
        }
        knn_rtree_time.push_back(knn_rtree_time_avg);
    }

    bool all_identical_results = true;
    if (knn_grid_results != knn_grid_fast_results &&
        knn_grid_results != knn_grid_compress_results)
        // knn_grid_results != knn_grid_disk_results)
        all_identical_results = false;
    for (int i = 0; i < index_n.size() - 1; i++)
        if (knn_grid_results[i] != knn_grid_results[i + 1])
            all_identical_results = false;
    if (knn_linear_scan_results != knn_grid_results[0])
        all_identical_results = false;
    if (knn_rtree_results != knn_grid_results[0])
        all_identical_results = false;
    if (all_identical_results)
        cout << "All methods gave identical results." << endl;
    else cout << "Results are different." << endl;

    if (strcmp(eval_data_path, "_") == 0)
        return 0;
    fstream eval_data_file;
    eval_data_file.open(eval_data_path, ios::out);
    if (!eval_data_file.is_open()) {
        cout << endl << "FAILED TO OPEN " << eval_data_path << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    eval_data_file << fixed;
    eval_data_file << "n,index loading time (normal)" << endl;
    for (int i = 0; i < grid_index_loading_time.size(); i++)
        eval_data_file << index_n[i] << "," << grid_index_loading_time[i] << endl;
    eval_data_file << "n,index loading time (compress)" << endl;
    for (int i = 0; i < grid_index_compress_loading_time.size(); i++)
        eval_data_file << index_n[i] << "," << grid_index_compress_loading_time[i] << endl;
    // eval_data_file << "n,index loading time (disk)" << endl;
    // for (int i = 0; i < grid_index_disk_loading_time.size(); i++)
    //     eval_data_file << index_n[i] << "," << grid_index_disk_loading_time[i] << endl;
    eval_data_file << "M,r-tree loading time" << endl << RTREE_MAX_ELEMENTS << "," << rtree_loading_time << endl;
    eval_data_file << "k,knn_linear_scan";
    for (int i = 0; i < index_n.size(); i++) {
        eval_data_file << ",knn_grid (n=" << index_n[i] << ")";
        eval_data_file << ",knn_grid_fast (n=" << index_n[i] << ")";
        eval_data_file << ",knn_grid_compress (n=" << index_n[i] << ")";
        // eval_data_file << ",knn_grid_disk (n=" << index_n[i] << ")";
    }
    eval_data_file << ",knn_rtree (M=" << RTREE_MAX_ELEMENTS << ")" << endl;
    for (int k = 0; k < K_MAX; k++) {
        eval_data_file << k + 1 << "," << knn_linear_scan_time[k];
        for (int i = 0; i < index_n.size(); i++) {
            eval_data_file << "," << knn_grid_time[i][k];
            eval_data_file << "," << knn_grid_fast_time[i][k];
            eval_data_file << "," << knn_grid_compress_time[i][k];
            // eval_data_file << "," << knn_grid_disk_time[i][k];
        }
        eval_data_file << "," << knn_rtree_time[k] << endl;
    }
    eval_data_file.close();
    return 0;
}
*/

int main(int argc, char** argv) {
    if (argc != 7) {
        cout << "Usage: " << argv[0] << " X Y DATA_PATH_NEW INDEX_PATH K N" << endl;
    
        // X(double): the latitude of the query point q
        // Y(double): the longitude of the query point q
        // DATA_PATH_NEW(char *): the file path of dataset you generated without duplicates
        // INDEX_PATH(char *): the file path of the grid index
        // K(integer): the k value for k-NN search
        // N(integer): the grid index size
    
        return -1;
    }
    clock_t s, t;
    cout << "Loading the preprocessed dataset from " << argv[3] << endl;
    s = clock();
    list<SpatialPoint> spatial_data = preload_dataset(argv[3]);
    t = clock();
    cout << "   Loading time: " << double(t - s) << endl;
    
    cout << "\nLoading the grid index from " << argv[4] << " (normal)" << endl;
    s = clock();
    GridIndex grid_index(argv[4], atoi(argv[6]));
    t = clock();
    cout << "   Loading time: " << double(t - s) << endl;
    
    cout << "\nLoading the grid index from " << argv[4] << " (compress)" << endl;
    s = clock();
    GridIndexCompress grid_index_compress(argv[4], atoi(argv[6]));
    t = clock();
    cout << "   Loading time: " << double(t - s) << endl;
    
    cout << "\nLoading the grid index from " << argv[4] << " (disk)" << endl;
    s = clock();
    GridIndexDisk grid_index_disk(argv[4], atoi(argv[6]), BUFFER_RATE);
    t = clock();
    cout << "   Loading time: " << double(t - s) << endl;
    
    cout << "\nLoading the preprocessed dataset from " << argv[3] << " building an R-tree" << endl;
    s = clock();
    RTree rtree = preload_rtree(argv[3]);
    t = clock();
    cout << "   Loading time: " << double(t - s) << endl;
    
    cout << endl;
    
    s = clock();
    string knn_linear_result = knn_linear_scan_preload(spatial_data , atof(argv[1]), atof(argv[2]), atoi(argv[5]));
    t = clock();
    cout << "Linear scan results:\n" << knn_linear_result << endl;
    cout << "Linear scan time: " << (double)(t - s) << endl;
    cout << endl;

    s = clock();
    string knn_grid_result = grid_index.knn_search(atof(argv[1]), atof(argv[2]), atoi(argv[5]));
    t = clock();
    cout << "Grid index search results:\n" << knn_grid_result << endl;
    cout << "Grid index search time: " << (double)(t - s) << endl;
    cout << endl;
    
    s = clock();
    string knn_grid_fast_result = grid_index.knn_search_fast(atof(argv[1]), atof(argv[2]), atoi(argv[5]));
    t = clock();
    cout << "Grid index (fast) search results:\n" << knn_grid_fast_result << endl;
    cout << "Grid index (fast) search time: " << (double)(t - s) << endl;
    cout << endl;
    
    s = clock();
    string knn_grid_compress_result = grid_index_compress.knn_search(atof(argv[1]), atof(argv[2]), atoi(argv[5]));
    t = clock();
    cout << "Grid index (compress) search results:\n" << knn_grid_compress_result << endl;
    cout << "Grid index (compress) search time: " << (double)(t - s) << endl;
    cout << endl;

    s = clock();
    string knn_rtree_result = knn_rtree_preload(rtree, atof(argv[1]), atof(argv[2]), atoi(argv[5]));
    t = clock();
    cout << "R-tree search results:\n" << knn_rtree_result << endl;
    cout << "R-tree search time: " << (double)(t - s) << endl;
    cout << endl;

    // sometimes segmentation faults happen with knn_grid_disk,
    // so it is the last one to execute.
    s = clock();
    string knn_grid_disk_result = grid_index_disk.knn_search(atof(argv[1]), atof(argv[2]), atoi(argv[5]));
    t = clock();
    cout << "Grid index (disk) search results:\n" << knn_grid_disk_result << endl;
    cout << "Grid index (disk) search time: " << (double)(t - s) << endl;
    cout << endl;

    if (knn_linear_result == knn_grid_result &&
        knn_grid_result == knn_grid_fast_result &&
        knn_grid_fast_result == knn_grid_compress_result &&
        knn_grid_compress_result == knn_grid_disk_result &&
        knn_grid_disk_result == knn_rtree_result)
        cout << "All search methods return the same results." << endl;
    else cout << "Results are different." << endl;

    return 0;
}
