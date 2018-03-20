

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include <cmath>
#include <memory>
#include <time.h>


struct Client {
    Client(int id, int demand, double x, double y)
        : id(id), demand(demand), x(x), y(y) { }

    int id;
    int demand;
    double x, y;
};


class Route {
public:
    Route(int capacity, std::shared_ptr<std::vector<std::vector<double>>> map)
        :   capacity_(capacity), numHop_(2), map_(map) {
        order_.push_back(0);
        order_.push_back(0);
    }

    bool canServe(int demand) {
        return (demand <= capacity_) ? true : false;
    }

    void appendClient(int id, int demand) {
        order_.insert(order_.begin() + numHop_ - 1, id);
        capacity_ -= demand;
        ++numHop_;
    }

    void insertClient(int pos, int id, int demand) {
        order_.insert(order_.begin() + pos, id);
        capacity_ -= demand;
        ++numHop_;
    }

    void eraseClient(int pos, int id, int demand) {
        order_.erase(order_.begin() + pos);
        capacity_ += demand;
        --numHop_;
    }

    double countCost() {
        double total = 0;
        for (int i = 1 ; i < numHop_ ; ++i) {
            int src = order_[i - 1];
            int dst = order_[i];
            total += (*map_)[src][dst];
        }
        return total;
    }

    int getHopCount() {
        return numHop_;
    }

    std::vector<int> getOrder() {
        return order_;
    }

private:
    int capacity_;
    int numHop_;

    std::shared_ptr<std::vector<std::vector<double>>> map_;
    std::vector<int> order_;
};


class Vrp {
public:
    Vrp(const std::string& fileName);

    ~Vrp() = default;

    Vrp(const Vrp&) = delete;
    Vrp(Vrp&&) = delete;
    Vrp& operator=(const Vrp&) = delete;
    Vrp& operator=(Vrp&&) = delete;

public:
    void run();

private:
    int numClient_;
    int numVehicle_;
    int capacity_;

    std::vector<Client> clients_;
    std::shared_ptr<std::vector<std::vector<double>>> map_;

    double curCost_;
    double optCost_;
    std::vector<Route> curOperation_;
    std::vector<Route> optOperation_;
};


Vrp::Vrp(const std::string& fileName) {

    std::ifstream input(fileName);
    input >> numClient_ >> numVehicle_ >> capacity_;

    for (int i = 0 ; i < numClient_ ; ++i) {
        int demand;
        double x, y;
        input >> demand >> x >> y;
        clients_.push_back(Client(i, demand, x, y));
    }

    map_ = std::make_shared<std::vector<std::vector<double>>>(numClient_);
    for (int i = 0 ; i < numClient_ ; ++i) {
        (*map_)[i].reserve(numClient_);
    }

    for (int i = 0 ; i < numClient_ ; ++i) {
        for (int j = 0 ; j < numClient_ ; ++j) {
            int srcX = clients_[i].x;
            int srcY = clients_[i].y;
            int dstX = clients_[j].x;
            int dstY = clients_[j].y;

            double distance = std::sqrt(
                std::pow(srcX - dstX, 2) + std::pow(srcY - dstY, 2));
            (*map_)[i][j] = (*map_)[j][i] = distance;
        }
    }

    for (int i = 0 ; i < numVehicle_ ; ++i) {
        curOperation_.push_back(Route(capacity_, map_));
    }

    // Prepare the initial client groups.
    int routeIdx = 0;
    int capacity = capacity_;

    for (int i = 1 ; i < numClient_ ; ++i) {
        int demand = clients_[i].demand;
        capacity -= demand;
        if (capacity > 0) {
            curOperation_[routeIdx].appendClient(i, demand);
        } else {
            ++routeIdx;
            curOperation_[routeIdx].appendClient(i, demand);
            capacity = capacity_ - demand;
        }
    }

    curCost_ = 0;
    for (int i = 0 ; i < numVehicle_ ; ++i) {
        curCost_ += curOperation_[i].countCost();
    }
    optCost_ = curCost_;

    optOperation_ = curOperation_;
}

void Vrp::run() {

    std::cout << optCost_ << " 0" << std::endl;
    for (int i = 0 ; i < numVehicle_ ; ++i) {
        const auto& order = optOperation_[i].getOrder();
        for (int id : order) {
            std::cout << id << ' ';
        }
        std::cout << std::endl;
    }
}


int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Please specify: "
                     "INPUT_PATH" << std::endl;
        return -1;
    }

    std::string fileName(argv[1]);
    Vrp solver(fileName);
    solver.run();

    return 0;
}
