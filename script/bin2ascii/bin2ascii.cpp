#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <filesystem>

namespace fs = std::filesystem;

typedef double real;

void load_tecplot3D_serial_binary_to_ascii(const char* binary_filename, const char* ascii_filename,
    int& nx, int& ny, int& nz, real dx, real dy, real dz)
{
    std::ifstream file(binary_filename, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening binary file for reading: " << binary_filename << std::endl;
        return;
    }

    file.read((char*)&nx, sizeof(int));
    file.read((char*)&ny, sizeof(int));
    file.read((char*)&nz, sizeof(int));

    int totalPoints = nx * ny * nz;
    std::vector<real> rhg(totalPoints);
    std::vector<real> uxg(totalPoints);
    std::vector<real> uyg(totalPoints);
    std::vector<real> uzg(totalPoints);

    real x, y, z;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int ijk = k * ny * nx + j * nx + i;
                file.read((char*)&x, sizeof(real));
                file.read((char*)&y, sizeof(real));
                file.read((char*)&z, sizeof(real));
                file.read((char*)&rhg[ijk], sizeof(real));
                file.read((char*)&uxg[ijk], sizeof(real));
                file.read((char*)&uyg[ijk], sizeof(real));
                file.read((char*)&uzg[ijk], sizeof(real));
            }
        }
    }
    file.close();

    std::ofstream ascii_file(ascii_filename);
    if (!ascii_file.is_open()) {
        std::cerr << "Error opening ASCII file for writing: " << ascii_filename << std::endl;
        return;
    }

    ascii_file << "TITLE = \"Example: Simple 3D-Field Data\"" << std::endl;
    ascii_file << "VARIABLES = \"X\", \"Y\", \"Z\", \"RHO\", \"UX\", \"UY\", \"UZ\"" << std::endl;
    ascii_file << "ZONE I=" << nx << ", J=" << ny << ", K=" << nz << " F=POINT" << std::endl;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int ijk = k * ny * nx + j * nx + i;
                ascii_file << (i + 0.5) * dx << " "
                           << (j + 0.5) * dy << " "
                           << (k + 0.5) * dz << " "
                           << rhg[ijk] << " "
                           << uxg[ijk] << " "
                           << uyg[ijk] << " "
                           << uzg[ijk] << std::endl;
            }
        }
    }
    ascii_file.close();
}

void convert_all_bin_to_tec(const std::string& input_folder, const std::string& output_folder, real dx, real dy, real dz)
{
    std::cout << "Input Folder: " << input_folder << std::endl;
    std::cout << "Output Folder: " << output_folder << std::endl;
    std::cout << "==================================================" << std::endl;

    if (!fs::exists(output_folder)) {
        fs::create_directories(output_folder);
    }

    for (const auto& entry : fs::directory_iterator(input_folder)) {
        if (entry.is_regular_file() && entry.path().extension() == ".bin") {
            std::string binary_filename = entry.path().string();
            std::string output_filename = (fs::path(output_folder) / entry.path().stem()).string() + ".tec";

            std::cout << "Processing: " << binary_filename << std::endl;
            std::cout << "Output to:  " << output_filename << std::endl;

            int nx, ny, nz;
            load_tecplot3D_serial_binary_to_ascii(binary_filename.c_str(), output_filename.c_str(), nx, ny, nz, dx, dy, dz);
        }
    }
}

int main() {
    std::string input_folder = "D:/FromTAIYI/DUGKS-PFM-GPU/100Channel/F3D/";
    std::string output_folder = "D:/FromTAIYI/DUGKS-PFM-GPU/100Channel/F3D/tec_output/";
    real dx = 0.482188, dy = 0.482188, dz = 0.482188;

    std::cout << "Grid spacing: dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;

    convert_all_bin_to_tec(input_folder, output_folder, dx, dy, dz);

    return 0;
}
