#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include <cmath>
#include <time.h>


class TabuList {
public:
    TabuList(int capacity)
        : size_(0), capacity_(capacity) { }
    ~TabuList() = default;

    TabuList(const TabuList&) = delete;
    TabuList(TabuList&&) = delete;
    TabuList& operator=(const TabuList&) = delete;
    TabuList& operator=(TabuList&&) = delete;

public:
    bool find(int key) {
        auto iter = set_.find(key);
        return iter == set_.cend() ? false : true;
    }

    void insert(int key) {
        if (size_ == capacity_) {
            auto item = queue_.top();
            queue_.pop();
            set_.erase(item);
            --size_;
        }

        queue_.push(key);
        set_.insert(key);
        ++size_;
    }

private:
    int size_;
    int capacity_;

    std::priority_queue<int> queue_;
    std::unordered_multiset<int> set_;
};


class TwoOpt {
public:
    TwoOpt(const std::string& fileName, double biasRate, int dropRate);
    ~TwoOpt() = default;

    TwoOpt(const TwoOpt&) = delete;
    TwoOpt(TwoOpt&&) = delete;
    TwoOpt& operator=(const TwoOpt&) = delete;
    TwoOpt& operator=(TwoOpt&&) = delete;

public:
    void run();
    std::pair<double, std::vector<int>> getAnswer() {
        return make_pair(optTourLen_, optTour_);
    }

private:
    void swapEdge(int k, int h);

private:
    int numNode_;
    int numTrial_;

    double biasRate_;
    int dropRate_;

    double ansTourLen_;
    double optTourLen_;
    double curTourLen_;

    std::vector<std::pair<double, double>> map_;
    std::vector<std::vector<double>> distance_;

    std::vector<int> ansTour_;
    std::vector<int> optTour_;
    std::vector<int> curTour_;
};


TwoOpt::TwoOpt(const std::string& fileName, double biasRate, int dropRate)
    : numTrial_(1), biasRate_(biasRate), dropRate_(dropRate) {

    std::ifstream input(fileName);
    input >> numNode_;

    map_.resize(numNode_);
    distance_.resize(numNode_);
    for (int i = 0 ; i < numNode_ ; ++i) {
        distance_[i].resize(numNode_);
        input >> map_[i].first >> map_[i].second;
    }

    for (int i = 0 ; i < numNode_ ; ++i) {
        for (int j = i + 1 ; j < numNode_ ; ++j) {
            double xSqr = pow(map_[j].first - map_[i].first, 2);
            double ySqr = pow(map_[j].second - map_[i].second, 2);
            distance_[i][j] = distance_[j][i] = sqrt(xSqr + ySqr);
        }
    }

    // Construct the initial tour
    optTourLen_ = 0;
    curTour_.resize(numNode_ + 1);
    optTour_.resize(numNode_ + 1);
    optTour_[0] = 0;

    for (int i = 1 ; i < numNode_ ; ++i) {
        optTour_[i] = i;
        optTourLen_ += distance_[i][i - 1];
    }

    optTour_[numNode_] = 0;
    optTourLen_ += distance_[0][numNode_ - 1];

    ansTourLen_ = optTourLen_;
}

void TwoOpt::run() {

    srand(time(NULL));

    int numFail = 0;
    TabuList tabu(numNode_ * 100);

    while (true) {
        // Continuously swap edges until no improvement can be found.
        for (int i = 1 ; i < numNode_ - 3 ; ++i) {
            for (int j = i + 1 ; j < numNode_ ; ++j) {

                swapEdge(i, j);
                if (tabu.find(curTourLen_)) {
                    continue;
                }

                if (curTourLen_ < optTourLen_) {
                    optTourLen_ = curTourLen_;
                    optTour_ = curTour_;

                    if (curTourLen_ < ansTourLen_) {
                        ansTourLen_ = curTourLen_;
                        ansTour_ = curTour_;
                    }

                    numFail = 0;
                    continue;
                }

                double diff = curTourLen_ - optTourLen_;
                if (diff / optTourLen_ < biasRate_) {
                    auto pick = rand() % dropRate_;
                    if (pick == 1) {
                        optTourLen_ = curTourLen_;
                        optTour_ = curTour_;
                        numFail = 0;
                    }
                    continue;
                }

                tabu.insert(curTourLen_);

                ++numFail;
                if (numFail == numNode_ * numNode_) {
                    return;
                }
            }
        }
    }
}

void TwoOpt::swapEdge(int k, int h) {

    for (int i = 0 ; i < k ; ++i) {
        curTour_[i] = optTour_[i];
    }

    for (int i = k, j = h ; i <= h ; ++i, --j) {
        curTour_[i] = optTour_[j];
    }

    for (int i = h + 1 ; i <= numNode_ ; ++i) {
        curTour_[i] = optTour_[i];
    }

    curTourLen_ = 0;
    for (int i = 1 ; i <= numNode_ ; ++i) {
        int src = curTour_[i - 1];
        int dst = curTour_[i];
        curTourLen_ += distance_[src][dst];
    }
}


int main(int argc, char** argv) {

    if (argc != 4) {
        std::cerr << "Please specify the input data path, bias rate, and random drop rate" << std::endl;
        return -1;
    }

    std::string fileName(argv[1]);
    double biasRate = atof(argv[2]);
    int dropRate = atoi(argv[3]);

    TwoOpt solver(fileName, biasRate, dropRate);
    solver.run();

    auto pair = solver.getAnswer();
    std::cout.precision(4);
    std::cout << pair.first << std::endl;
    for (int i = 0 ; i < pair.second.size() - 1; ++i) {
        std::cout << pair.second[i] << ' ';
    }
    std::cout << std::endl;

    return 0;
}