#include "happly.h"
#include <iostream>

int main() {
    // Construct a data object by reading from file
    happly::PLYData data("test_data.ply");
    auto elem_names = data.getElementNames();
    for (const auto& name : elem_names) {
        std::cout << name << std::endl;
    }
    auto positions = data.getVertexPositions();
    std::cout << positions.size() << std::endl;
    auto colors = data.getVertexColors();
}