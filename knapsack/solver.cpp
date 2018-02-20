#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <algorithm>


struct Item {
    Item(int index, int value, int weight)
        :   index(index),
            value(value),
            weight(weight),
            ratio(static_cast<double>(value) / weight) { }

    int index;
    int value;
    int weight;
    double ratio;
};


class BranchAndBound {
public:
    BranchAndBound(const std::string& fileName);
    ~BranchAndBound() = default;

    BranchAndBound(const BranchAndBound&) = delete;
    BranchAndBound(BranchAndBound&&) = delete;
    BranchAndBound& operator=(const BranchAndBound&) = delete;
    BranchAndBound& operator=(BranchAndBound&&) = delete;

public:
    void run(int depth);
    std::pair<int, std::vector<bool>> getAnswer();

private:
    int getEstimation(int depth);

private:
    int numItems_;
    int maxCapacity_;
    int curCapacity_;
    int optimal_;
    int local_;

    std::vector<bool> optimalMap_;
    std::vector<bool> localMap_;

    std::vector<std::shared_ptr<Item>> sortedItems_;
};


BranchAndBound::BranchAndBound(const std::string& fileName)
    :   optimal_(-1), local_(0) {

    std::ifstream input(fileName);
    input >> numItems_ >> maxCapacity_;
    curCapacity_ = maxCapacity_;

    localMap_.resize(numItems_);

    for (int i = 0 ; i < numItems_ ; ++i) {
        int value, weight;
        input >> value >> weight;
        sortedItems_.push_back(std::make_shared<Item>(i, value, weight));
    }

    auto compare = [] (std::shared_ptr<Item> lhs, std::shared_ptr<Item> rhs) {
        return lhs->ratio > rhs->ratio;
    };
    std::sort(sortedItems_.begin(), sortedItems_.end(), compare);
}

void BranchAndBound::run(int depth) {

    if (depth == numItems_) {
        return;
    }

    auto item = sortedItems_[depth];

    localMap_[depth] = true;
    local_ += item->value;
    curCapacity_ -= item->weight;

    if (curCapacity_ >= 0 && local_ > optimal_) {
        optimal_ = local_;
        optimalMap_ = localMap_;
    }

    int estimation = getEstimation(depth + 1);
    if (estimation >= optimal_) {
        run(depth + 1);
    }

    localMap_[depth] = false;
    local_ -= item->value;
    curCapacity_ += item->weight;

    estimation = getEstimation(depth + 1);
    if (estimation >= optimal_) {
        run(depth + 1);
    }
}

int BranchAndBound::getEstimation(int depth) {

    if (curCapacity_ < 0) {
        return 0;
    }

    int estimation = local_;
    int capacity = curCapacity_;

    while (depth < numItems_) {
        auto item = sortedItems_[depth];
        int weight = item->weight;
        if (weight > capacity) {
            break;
        }

        estimation += item->value;
        capacity -= weight;
        ++depth;
    }

    if (depth < numItems_) {
        auto item = sortedItems_[depth];
        estimation += item->value * capacity / item->weight;
    }

    return estimation;
}

std::pair<int, std::vector<bool>> BranchAndBound::getAnswer() {
    std::vector<bool> answerMap(numItems_, false);

    for (int i = 0 ; i < numItems_ ; ++i) {
        bool select = optimalMap_[i];
        if (!select) {
            continue;
        }

        int index = sortedItems_[i]->index;
        answerMap[index] = true;
    }

    return std::make_pair(optimal_, answerMap);
}


int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Please specify the input data path" << std::endl;
        return -1;
    }

    std::string fileName(argv[1]);
    BranchAndBound solver(fileName);
    solver.run(0);

    auto answer = solver.getAnswer();
    std::cout << answer.first << " 0" << std::endl;
    for (auto bit : answer.second) {
        std::cout << bit << ' ';
    }
    std::cout << std::endl;

    return 0;
}