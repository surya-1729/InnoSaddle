#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

void convert_obj_to_off(std::string input_obj_file, std::string output_off_file) {
    std::ifstream input(input_obj_file);
    if (!input.is_open()) {
        std::cerr << "Could not open the file: " << input_obj_file << std::endl;
        return;
    }

    std::vector<std::string> vertices;
    std::vector<std::vector<int>> faces;

    std::string line;
    while (std::getline(input, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v") {
            vertices.push_back(line.substr(2));
        } else if (prefix == "f") {
            std::vector<int> face;
            int index;
            while (iss >> index) {
                face.push_back(index);
            }
            faces.push_back(face);
        }
    }
    input.close();

    std::ofstream output(output_off_file);
    if (!output.is_open()) {
        std::cerr << "Could not open the file: " << output_off_file << std::endl;
        return;
    }

    output << "OFF\n";
    output << vertices.size() << " " << faces.size() << " 0\n";

    for (const auto& vertex : vertices) {
        output << vertex << "\n";
    }

    for (const auto& face : faces) {
        output << face.size() << " ";
        for (const auto& index : face) {
            output << (index - 1) << " ";
        }
        output << "\n";
    }
    output.close();
}

int main() {
    convert_obj_to_off("./data/processed_data/0004_piri.obj", "./data/processed_data/0004_piri.off");
    return 0;
}
