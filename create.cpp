#include <fstream>
#include <string>

void strip_unsupported_obj_data(const std::string& input_path, const std::string& output_path) {
    std::ifstream input_file(input_path);
    std::ofstream output_file(output_path);

    std::string line;
    while (std::getline(input_file, line)) {
        if (line.substr(0, 2) == "v " || line.substr(0, 2) == "f ") {
            output_file << line << '\n';
        }
    }
}


int main() {
    std::string input_path = "./data/processed_data/0004_piri.obj";
    std::string output_path = "./data/processed_data/0004_piri_stripped.obj";
    strip_unsupported_obj_data(input_path, output_path);

    // continue with your code...
    return 0;
}
