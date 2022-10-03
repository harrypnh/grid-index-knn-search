#ifndef GRIDINDEX_H
#define GRIDINDEX_H

#include <vector>
#include <list>
#include <string>
#include <unordered_map>

#define NOT_IN_MAXBOX -1
#define CELL_NOT_FOUND -2
#define PRECISION 10
#define EPSILON 1e-10

using namespace std;

bool is_equal(double x, double y);

class SpatialPoint {
    private:
        double latitude;
        double longitude;
        int location_id;
        double query_distance;
        int cell_i, cell_j;
        bool located_on_grid;
    public:
        SpatialPoint(double latitude, double longitude, int location_id);
        SpatialPoint(double latitude, double longitude, int location_id, double query_distance);
        SpatialPoint(const SpatialPoint &p);
        void operator=(const SpatialPoint &p);
        bool is_located_on_grid();
        double get_latitude();
        double get_longitude();
        int get_location_id();
        int get_cell_i();
        int get_cell_j();
        double get_query_distance();
        void set_latitude(double latitude);
        void set_longitude(double longitude);
        void set_location_id(int location_id);
        void set_cell_container(int cell_i, int cell_j);
        void calculate_query_distance(SpatialPoint q);
};

bool compare_spatial_point(SpatialPoint p1, SpatialPoint p2);

bool compare_query_distance(SpatialPoint p1, SpatialPoint p2);

class GridCell {
    private:
        list<SpatialPoint> cell_data;
        int cell_i, cell_j;
        double x_min, x_max, y_min, y_max;
        double query_distance;
    public:
        GridCell() = default;
        GridCell(int cell_i, int cell_j, double x_min, double x_max, double y_min, double y_max);
        GridCell(const GridCell &c);
        void operator=(const GridCell &c);
        list<SpatialPoint> get_cell_data();
        int get_cell_i();
        int get_cell_j();
        double get_query_distance();
        void set_query_distance(SpatialPoint q);
        void add_spatial_point(SpatialPoint p);
};

class GridIndexBase {
    protected:
        int n;
        double cell_range_x;
        double cell_range_y;
        double latitude_min;
        double latitude_max;
        double longitude_min;
        double longitude_max;

    public:
        void locate_point(SpatialPoint &p);
        virtual GridCell get_cell(int cell_i, int cell_j) = 0;
        int key_function(int cell_i, int cell_j);
        string knn_search(double x, double y, int k);
        string knn_search_fast(double x, double y, int k);
};

class GridIndex : public GridIndexBase {
    private:
        vector<vector<GridCell>> grid_index;
        
    public:
        GridIndex(char* index_path, int n);
        GridCell get_cell(int cell_i, int cell_j);
};

class GridIndexCompress : public GridIndexBase {
    private:
        unordered_map<int, GridCell> grid_index_compress;
    
    public:
        GridIndexCompress(char* index_path, int n);
        GridCell get_cell(int cell_i, int cell_j);
};

class GridIndexDisk : public GridIndexBase {
    private:
        int buffer_size;
        char* index_path;
        list<GridCell> buffer;
        unordered_map<int, list<GridCell>::iterator> hash_table;

    public:
        GridIndexDisk(char* index_path, int n, double buffer_rate);
        GridCell get_cell(int cell_i, int cell_j);
};

#endif