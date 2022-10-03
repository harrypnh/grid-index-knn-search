#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <algorithm>
#include <queue>
#include "GridIndex.h"

using namespace std;

bool is_equal(double x, double y) { return fabs(x - y) < EPSILON; }

SpatialPoint::SpatialPoint(double latitude, double longitude, int location_id) {
    this->latitude = latitude;
    this->longitude = longitude;
    this->location_id = location_id;
    query_distance = DBL_MAX;
    located_on_grid = false;
}

SpatialPoint::SpatialPoint(double latitude, double longitude, int location_id, double query_distance) {
    this->latitude = latitude;
    this->longitude = longitude;
    this->location_id = location_id;
    this->query_distance = query_distance;
    located_on_grid = false;
}

SpatialPoint::SpatialPoint(const SpatialPoint &p) {
    this->latitude = p.latitude;
    this->longitude = p.longitude;
    this->location_id = p.location_id;
    this->query_distance = p.query_distance;
    this->cell_i = p.cell_i;
    this->cell_j = p.cell_j;
    this->located_on_grid = p.located_on_grid;
}

void SpatialPoint::operator=(const SpatialPoint &p) {
    this->latitude = p.latitude;
    this->longitude = p.longitude;
    this->location_id = p.location_id;
    this->query_distance = p.query_distance;
    this->cell_i = p.cell_i;
    this->cell_j = p.cell_j;
    this->located_on_grid = p.located_on_grid;
}

bool SpatialPoint::is_located_on_grid() { return located_on_grid; }

double SpatialPoint::get_latitude() { return latitude; }

double SpatialPoint::get_longitude() { return longitude; }

int SpatialPoint::get_location_id() { return location_id; }

double SpatialPoint::get_query_distance() { return query_distance; }

int SpatialPoint::get_cell_i() { return cell_i; }

int SpatialPoint::get_cell_j() { return cell_j; }

void SpatialPoint::set_latitude(double latitude) { this->latitude = latitude; }

void SpatialPoint::set_longitude(double longitude) { this->longitude = longitude; }

void SpatialPoint::set_location_id(int location_id) { this->location_id = location_id; }

void SpatialPoint::set_cell_container(int cell_i, int cell_j) {
    this->cell_i = cell_i;
    this->cell_j = cell_j;
    if (cell_i == NOT_IN_MAXBOX || cell_j == NOT_IN_MAXBOX)
        located_on_grid = false;
    else
        located_on_grid = true;
}

void SpatialPoint::calculate_query_distance(SpatialPoint q) {
    query_distance = sqrt(pow(q.get_latitude() - latitude, 2.0) + pow(q.get_longitude() - longitude, 2.0));
}

bool compare_spatial_point(SpatialPoint p1, SpatialPoint p2) {
    if (p1.get_latitude() == p2.get_latitude() && p1.get_longitude() == p2.get_longitude())
        return p1.get_location_id() < p2.get_location_id();
    if (p1.get_latitude() == p2.get_latitude())
        return p1.get_longitude() < p2.get_longitude();
    return p1.get_latitude() < p2.get_latitude();
}

bool compare_query_distance(SpatialPoint p1, SpatialPoint p2) {
    return p1.get_query_distance() < p2.get_query_distance();
}

GridCell::GridCell(int cell_i, int cell_j, double x_min, double x_max, double y_min, double y_max) {
    this->cell_i = cell_i;
    this->cell_j = cell_j;
    this->x_min = x_min;
    this->x_max = x_max;
    this->y_min = y_min;
    this->y_max = y_max;
}

GridCell::GridCell(const GridCell &c) {
    this->cell_data = c.cell_data;
    this->cell_i = c.cell_i;
    this->cell_j = c.cell_j;
    this->x_min = c.x_min;
    this->x_max = c.x_max;
    this->y_min = c.y_min;
    this->y_max = c.y_max;
    this->query_distance = c.query_distance;
}

void GridCell::operator=(const GridCell &c) {
    this->cell_data = c.cell_data;
    this->cell_i = c.cell_i;
    this->cell_j = c.cell_j;
    this->x_min = c.x_min;
    this->x_max = c.x_max;
    this->y_min = c.y_min;
    this->y_max = c.y_max;
    this->query_distance = c.query_distance;
}

list<SpatialPoint> GridCell::get_cell_data() { return cell_data; }

int GridCell::get_cell_i() { return cell_i; }

int GridCell::get_cell_j() { return cell_j; }

double GridCell::get_query_distance() { return query_distance; }

void GridCell::set_query_distance(SpatialPoint q) {
    if (!q.is_located_on_grid())
        query_distance = -1.0;
    else if (q.get_cell_i() == cell_i && q.get_cell_j() == cell_j)
        query_distance = 0.0;
    else {
        double x_q = q.get_latitude(), y_q = q.get_longitude();
        if (q.get_cell_i() == cell_i && q.get_cell_j() < cell_j)
            query_distance = x_min - x_q;
        else if (q.get_cell_i() == cell_i && q.get_cell_j() > cell_j)
            query_distance = x_q - x_max;
        else if (q.get_cell_j() == cell_j && q.get_cell_i() < cell_i)
            query_distance = y_min - y_q;
        else if (q.get_cell_j() == cell_j && q.get_cell_i() > cell_i)
            query_distance = y_q - y_max;
        else {
            vector<double> distance_to_corners;
            distance_to_corners.reserve(4);
            distance_to_corners.push_back(sqrt(pow(x_q - x_min, 2.0) + pow(y_q - y_min, 2.0)));
            distance_to_corners.push_back(sqrt(pow(x_q - x_max, 2.0) + pow(y_q - y_min, 2.0)));
            distance_to_corners.push_back(sqrt(pow(x_q - x_min, 2.0) + pow(y_q - y_max, 2.0)));
            distance_to_corners.push_back(sqrt(pow(x_q - x_max, 2.0) + pow(y_q - y_max, 2.0)));
            double dlow = distance_to_corners[0];
            for (int i = 1; i < distance_to_corners.size(); i++)
                if (distance_to_corners[i] < dlow)
                    dlow = distance_to_corners[i];
            query_distance = dlow;
        }
    }
}

void GridCell::add_spatial_point(SpatialPoint p) { cell_data.push_back(p); }

void GridIndexBase::locate_point(SpatialPoint &p) {
    int cell_i, cell_j;
    if (p.get_latitude() < latitude_min || p.get_latitude() > latitude_max ||
        p.get_longitude() < longitude_min || p.get_longitude() > longitude_max) {
        p.set_cell_container(NOT_IN_MAXBOX, NOT_IN_MAXBOX);
        return;
    }
    double i_approx = (p.get_longitude() - longitude_min) / cell_range_y;
    double j_approx = (p.get_latitude() - latitude_min) / cell_range_x;

    double i_approx_floor = floor(i_approx), j_approx_floor = floor(j_approx);
    if (is_equal(0.0, i_approx_floor))
        cell_i = 0;
    else if (is_equal(i_approx, i_approx_floor))
        cell_i = int(i_approx_floor) - 1;
    else if (i_approx > i_approx_floor)
        cell_i = int(i_approx_floor);

    if (is_equal(0.0, j_approx_floor))
        cell_j = 0;
    else if (is_equal(j_approx, j_approx_floor))
        cell_j = int(j_approx_floor) - 1;
    else if (j_approx > j_approx_floor)
        cell_j = int(j_approx_floor);
    p.set_cell_container(cell_i, cell_j);
    return;
}

int GridIndexBase::key_function(int cell_i, int cell_j) { return n * cell_i + cell_j; }

string GridIndexBase::knn_search(double x, double y, int k) {
    SpatialPoint query_point(x, y, 0, 0);
    locate_point(query_point);
    if (query_point.get_cell_i() == NOT_IN_MAXBOX || query_point.get_cell_j() == NOT_IN_MAXBOX) {
        cout << "QUERY POINT OUT OF MAXBOX." << endl;
        exit(-1);
    }
    list<SpatialPoint> knn;
    double knn_max_distance = DBL_MAX;
    GridCell first_cell = get_cell(query_point.get_cell_i(), query_point.get_cell_j());
    list<SpatialPoint> cell_data = first_cell.get_cell_data();
    for (auto point_it = cell_data.begin(); point_it != cell_data.end(); point_it++) {
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
    int n_layers = max(max(n - 1 - query_point.get_cell_i(), query_point.get_cell_i()), max(n - 1 - query_point.get_cell_j(), query_point.get_cell_j()));
    for (int layer = 1; layer <= n_layers; layer++) {
        bool prune_entire_layer = true;
        list<GridCell> cell_list;
        if (query_point.get_cell_i() - layer >= 0) {
            for (int j = query_point.get_cell_j() - layer; j <= query_point.get_cell_j() + layer; j++) {
                if (j >= 0 && j < n) {
                    GridCell cell = get_cell(query_point.get_cell_i() - layer, j);
                    cell_list.push_back(cell);
                }
            }
        }
        if (query_point.get_cell_i() + layer < n) {
            for (int j = query_point.get_cell_j() - layer; j <= query_point.get_cell_j() + layer; j++) {
                if (j >= 0 && j < n) {
                    GridCell cell = get_cell(query_point.get_cell_i() + layer, j);
                    cell_list.push_back(cell);
                }
            }
        }
        if (query_point.get_cell_j() - layer >= 0) {
            for (int i = query_point.get_cell_i() - layer + 1; i <= query_point.get_cell_i() + layer - 1; i++) {
                if (i >= 0 && i < n) {
                    GridCell cell = get_cell(i, query_point.get_cell_j() - layer);
                    cell_list.push_back(cell);
                }
            }
        }
        if (query_point.get_cell_j() + layer < n) {
            for (int i = query_point.get_cell_i() - layer + 1; i <= query_point.get_cell_i() + layer - 1; i++) {
                if (i >= 0 && i < n) {
                    GridCell cell = get_cell(i, query_point.get_cell_j() + layer);
                    cell_list.push_back(cell);
                }
            }
        }
        for (auto cell_it = cell_list.begin(); cell_it != cell_list.end(); cell_it++) {
            cell_it->set_query_distance(query_point);
            double dlow = cell_it->get_query_distance();
            if (dlow > knn_max_distance && knn.size() == k)
                continue;
            prune_entire_layer = false;
            list<SpatialPoint> cell_data = cell_it->get_cell_data();
            for (auto point_it = cell_data.begin(); point_it != cell_data.end(); point_it++) {
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
        }
        if (prune_entire_layer && knn.size() == k)
            break;
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

string GridIndexBase::knn_search_fast(double x, double y, int k) {
    SpatialPoint query_point(x, y, 0, 0);
    locate_point(query_point);
    if (query_point.get_cell_i() == NOT_IN_MAXBOX || query_point.get_cell_j() == NOT_IN_MAXBOX) {
        cout << "QUERY POINT OUT OF MAXBOX." << endl;
        exit(-1);
    }
    vector<SpatialPoint> knn_container;
    knn_container.reserve(k + 1);
    priority_queue<SpatialPoint, vector<SpatialPoint>, decltype(&compare_query_distance)> knn(compare_query_distance, knn_container);
    GridCell first_cell = get_cell(query_point.get_cell_i(), query_point.get_cell_j());
    list<SpatialPoint> cell_data = first_cell.get_cell_data();
    for (auto point_it = cell_data.begin(); point_it != cell_data.end(); point_it++) {
        point_it->calculate_query_distance(query_point);
        knn.push(SpatialPoint(point_it->get_latitude(), point_it->get_longitude(), point_it->get_location_id(), point_it->get_query_distance()));
        if (knn.size() == k + 1)
            knn.pop();
    } 
    double knn_max_distance = DBL_MAX;
    if (!knn.empty()) {
        SpatialPoint farthest_nn = knn.top();
        knn_max_distance = farthest_nn.get_query_distance();
    }
    int n_layers = max(max(n - 1 - query_point.get_cell_i(), query_point.get_cell_i()), max(n - 1 - query_point.get_cell_j(), query_point.get_cell_j()));
    for (int layer = 1; layer <= n_layers; layer++) {
        bool prune_entire_layer = true;
        list<GridCell> cell_list;
        if (query_point.get_cell_i() - layer >= 0) {
            for (int j = query_point.get_cell_j() - layer; j <= query_point.get_cell_j() + layer; j++) {
                if (j >= 0 && j < n) {
                    GridCell cell = get_cell(query_point.get_cell_i() - layer, j);
                    cell_list.push_back(cell);
                }
            }
        }
        if (query_point.get_cell_i() + layer < n) {
            for (int j = query_point.get_cell_j() - layer; j <= query_point.get_cell_j() + layer; j++) {
                if (j >= 0 && j < n) {
                    GridCell cell = get_cell(query_point.get_cell_i() + layer, j);
                    cell_list.push_back(cell);
                }
            }
        }
        if (query_point.get_cell_j() - layer >= 0) {
            for (int i = query_point.get_cell_i() - layer + 1; i <= query_point.get_cell_i() + layer - 1; i++) {
                if (i >= 0 && i < n) {
                    GridCell cell = get_cell(i, query_point.get_cell_j() - layer);
                    cell_list.push_back(cell);
                }
            }
        }
        if (query_point.get_cell_j() + layer < n) {
            for (int i = query_point.get_cell_i() - layer + 1; i <= query_point.get_cell_i() + layer - 1; i++) {
                if (i >= 0 && i < n) {
                    GridCell cell = get_cell(i, query_point.get_cell_j() + layer);
                    cell_list.push_back(cell);
                }
            }
        }
        for (auto cell_it = cell_list.begin(); cell_it != cell_list.end(); cell_it++) {
            cell_it->set_query_distance(query_point);
            double dlow = cell_it->get_query_distance();
            if (dlow > knn_max_distance && knn.size() == k)
                continue;
            prune_entire_layer = false;
            list<SpatialPoint> cell_data = cell_it->get_cell_data();
            for (auto point_it = cell_data.begin(); point_it != cell_data.end(); point_it++) {
                point_it->calculate_query_distance(query_point);
                if (point_it->get_query_distance() < knn_max_distance && knn.size() == k)
                    knn.pop();
                if (knn.size() < k) {
                    knn.push(SpatialPoint(point_it->get_latitude(), point_it->get_longitude(), point_it->get_location_id(), point_it->get_query_distance()));
                    SpatialPoint farthest_nn = knn.top();
                    knn_max_distance = farthest_nn.get_query_distance();
                }
            }
        }
        if (prune_entire_layer && knn.size() == k)
            break;
    }
    list<SpatialPoint> results;
    while(!knn.empty()) {
        results.push_back(knn.top());
        knn.pop();
    }
    ostringstream results_ss;
    for (auto it = results.rbegin(); it != results.rend(); it++) {
        if (it != prev(results.rend()))
            results_ss << it->get_location_id() << ", ";
        else
            results_ss << it->get_location_id();
    }
    return results_ss.str();
}

GridIndex::GridIndex(char* index_path, int n) {
    this->n = n;
    fstream index_file;
    index_file.open(index_path, ios::in);
    if (!index_file.is_open()) {
        cout << "FAILED TO OPEN " << index_path << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    // assume every index file starts with a line containing x_min x_max y_min y_max
    string max_box;
    getline(index_file, max_box);
    istringstream(max_box) >> latitude_min >> latitude_max >> longitude_min >> longitude_max;
    cell_range_x = (latitude_max - latitude_min) / n;
    cell_range_y = (longitude_max - longitude_min) / n;
    // assume cells are in the ascending order of rows then columns.
    string raw_cell_data;
    grid_index.reserve(n);
    for (int i = 0; i < n; i++) {
        vector<GridCell> row_data;
        row_data.reserve(n);
        for (int j = 0; j < n; j++) {
            getline(index_file, raw_cell_data);
            istringstream raw_cd_sstream(raw_cell_data);
            string temp;
            raw_cd_sstream >> temp;
            raw_cd_sstream >> temp;
            double x_min = latitude_min + j * cell_range_x;
            double x_max = latitude_min + (j + 1) * cell_range_x;
            double y_min = longitude_min + i * cell_range_y;
            double y_max = longitude_min + (i + 1) * cell_range_y;
            GridCell cell(i, j, x_min, x_max, y_min, y_max);
            while (!raw_cd_sstream.eof()) {
                string raw_spatial_point;
                raw_cd_sstream >> raw_spatial_point;
                int location_id;
                double latitude, longitude;
                istringstream raw_sp_sstream(raw_spatial_point);
                getline(raw_sp_sstream, temp, '_');
                istringstream(temp) >> location_id;
                getline(raw_sp_sstream, temp, '_');
                istringstream(temp) >> latitude;
                getline(raw_sp_sstream, temp, '_');
                istringstream(temp) >> longitude;
                SpatialPoint point(latitude, longitude, location_id);
                point.set_cell_container(i, j);
                cell.add_spatial_point(point);
            }
            row_data.push_back(cell);
        }
        grid_index.push_back(row_data);
    }
}

GridCell GridIndex::get_cell(int cell_i, int cell_j) { return grid_index[cell_i][cell_j]; }

GridIndexCompress::GridIndexCompress(char* index_path, int n) {
    this->n = n;
    fstream index_file;
    index_file.open(index_path, ios::in);
    if (!index_file.is_open()) {
        cout << "FAILED TO OPEN " << index_path << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    // assume every index file starts with a line containing x_min x_max y_min y_max
    string max_box;
    getline(index_file, max_box);
    istringstream(max_box) >> latitude_min >> latitude_max >> longitude_min >> longitude_max;
    cell_range_x = (latitude_max - latitude_min) / n;
    cell_range_y = (longitude_max - longitude_min) / n;
    // assume cells are in the ascending order of rows then columns.
    string raw_cell_data;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            getline(index_file, raw_cell_data);
            istringstream raw_cd_sstream(raw_cell_data);
            string temp;
            raw_cd_sstream >> temp;
            raw_cd_sstream >> temp;
            if (raw_cd_sstream.eof())
                continue;
            double x_min = latitude_min + j * cell_range_x;
            double x_max = latitude_min + (j + 1) * cell_range_x;
            double y_min = longitude_min + i * cell_range_y;
            double y_max = longitude_min + (i + 1) * cell_range_y;
            GridCell cell(i, j, x_min, x_max, y_min, y_max);
            while (!raw_cd_sstream.eof()) {
                string raw_spatial_point;
                raw_cd_sstream >> raw_spatial_point;
                int location_id;
                double latitude, longitude;
                istringstream raw_sp_sstream(raw_spatial_point);
                getline(raw_sp_sstream, temp, '_');
                istringstream(temp) >> location_id;
                getline(raw_sp_sstream, temp, '_');
                istringstream(temp) >> latitude;
                getline(raw_sp_sstream, temp, '_');
                istringstream(temp) >> longitude;
                SpatialPoint point(latitude, longitude, location_id);
                point.set_cell_container(i, j);
                cell.add_spatial_point(point);
            }
            grid_index_compress.insert(make_pair(key_function(i, j), cell));
        }
    }
}

GridCell GridIndexCompress::get_cell(int cell_i, int cell_j) {
    int key = key_function(cell_i, cell_j);
    if (grid_index_compress.find(key) != grid_index_compress.end())
        return grid_index_compress[key];
    double x_min = latitude_min + cell_j * cell_range_x;
    double x_max = latitude_min + (cell_j + 1) * cell_range_x;
    double y_min = longitude_min + cell_i * cell_range_y;
    double y_max = longitude_min + (cell_i + 1) * cell_range_y;
    GridCell empty_cell(cell_i, cell_j, x_min, x_max, y_min, y_max);
    return empty_cell;
}

GridIndexDisk::GridIndexDisk(char* index_path, int n, double buffer_rate) {
    this->index_path = index_path;
    this->n = n;
    buffer_size = floor(buffer_rate * n * n);
    fstream index_file;
    index_file.open(index_path, ios::in);
    if (!index_file.is_open()) {
        cout << "FAILED TO OPEN " << index_path << " OR FILE NOT FOUND." << endl;
        exit(-1);
    }
    // assume every index file starts with a line containing x_min x_max y_min y_max
    string max_box;
    getline(index_file, max_box);
    istringstream(max_box) >> latitude_min >> latitude_max >> longitude_min >> longitude_max;
    cell_range_x = (latitude_max - latitude_min) / n;
    cell_range_y = (longitude_max - longitude_min) / n;
    string raw_cell_data;
    while(getline(index_file, raw_cell_data) && buffer.size() < buffer_size) {
        istringstream raw_cd_sstream(raw_cell_data);
        string temp;
        raw_cd_sstream >> temp;
        raw_cd_sstream >> temp;
        istringstream cell_number(temp);
        string i_str, j_str;
        int i, j;
        getline(cell_number, i_str, ',');
        getline(cell_number, j_str, ':');
        istringstream(i_str) >> i;
        istringstream(j_str) >> j;
        double x_min = latitude_min + j * cell_range_x;
        double x_max = latitude_min + (j + 1) * cell_range_x;
        double y_min = longitude_min + i * cell_range_y;
        double y_max = longitude_min + (i + 1) * cell_range_y;
        GridCell cell(i, j, x_min, x_max, y_min, y_max);
        while (!raw_cd_sstream.eof()) {
            string raw_spatial_point;
            raw_cd_sstream >> raw_spatial_point;
            int location_id;
            double latitude, longitude;
            istringstream raw_sp_sstream(raw_spatial_point);
            getline(raw_sp_sstream, temp, '_');
            istringstream(temp) >> location_id;
            getline(raw_sp_sstream, temp, '_');
            istringstream(temp) >> latitude;
            getline(raw_sp_sstream, temp, '_');
            istringstream(temp) >> longitude;
            SpatialPoint point(latitude, longitude, location_id);
            point.set_cell_container(i, j);
            cell.add_spatial_point(point);
        }
        buffer.push_front(cell);
        hash_table.insert(make_pair(key_function(i, j), buffer.begin()));
    }
    index_file.close();
}

GridCell GridIndexDisk::get_cell(int cell_i, int cell_j) {
    int key = key_function(cell_i, cell_j);
    GridCell cell;
    if (hash_table.find(key) == hash_table.end()) {
        fstream index_file;
        index_file.open(index_path, ios::in);
        if (!index_file.is_open()) {
            cout << "FAILED TO OPEN " << index_path << " OR FILE NOT FOUND." << endl;
            exit(-1);
        }
        int lines_to_skip = 1 + key;
        while (!index_file.eof() && lines_to_skip > 0) {
            index_file.ignore(numeric_limits<streamsize>::max(), '\n');
            lines_to_skip--;
        }
        if (index_file.eof() && lines_to_skip > 0) {
            cout << "INDEX FILE " << index_path << " IS CORRUPTED" << endl;
            exit(-1);
        }
        string raw_cell_data;
        getline(index_file, raw_cell_data);
        istringstream raw_cd_sstream(raw_cell_data);
        string temp;
        raw_cd_sstream >> temp;
        raw_cd_sstream >> temp;
        double x_min = latitude_min + cell_j * cell_range_x;
        double x_max = latitude_min + (cell_j + 1) * cell_range_x;
        double y_min = longitude_min + cell_i * cell_range_y;
        double y_max = longitude_min + (cell_i + 1) * cell_range_y;
        cell = GridCell(cell_i, cell_j, x_min, x_max, y_min, y_max);
        while (!raw_cd_sstream.eof()) {
            string raw_spatial_point;
            raw_cd_sstream >> raw_spatial_point;
            int location_id;
            double latitude, longitude;
            istringstream raw_sp_sstream(raw_spatial_point);
            getline(raw_sp_sstream, temp, '_');
            istringstream(temp) >> location_id;
            getline(raw_sp_sstream, temp, '_');
            istringstream(temp) >> latitude;
            getline(raw_sp_sstream, temp, '_');
            istringstream(temp) >> longitude;
            SpatialPoint point(latitude, longitude, location_id);
            point.set_cell_container(cell_i, cell_j);
            cell.add_spatial_point(point);
        }
        index_file.close();
        if (buffer.size() == buffer_size) {
            hash_table.erase(key_function(buffer.back().get_cell_i(), buffer.back().get_cell_j()));
            buffer.pop_back();
        }
    }
    else {
        cell = *(hash_table[key]);
        buffer.erase(hash_table[key]);
        hash_table.erase(key);
    }
    buffer.push_front(cell);
    hash_table.insert(make_pair(key, buffer.begin()));
    return cell;
}
