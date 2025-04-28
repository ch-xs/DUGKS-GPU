#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <filesystem>  // C++17 filesystem library

namespace fs = std::filesystem;

// Assume 'real' is of type float
typedef double real;

void load_tecplot3D_serial_binary_to_ascii(const char* binary_filename, const char* ascii_filename,
    int& nx, int& ny, int& nz, real dx, real dy, real dz)
{
    std::ifstream file(binary_filename, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening binary file for reading: " << binary_filename << std::endl;
        return;
    }

    // Read grid size
    file.read((char*)&nx, sizeof(int));
    file.read((char*)&ny, sizeof(int));
    file.read((char*)&nz, sizeof(int));

    // Calculate total data points
    int totalPoints = nx * ny * nz;

    // Allocate space
    std::vector<real> rhg(totalPoints);
    std::vector<real> uxg(totalPoints);
    std::vector<real> uyg(totalPoints);
    std::vector<real> uzg(totalPoints);

    // Read data (7 real numbers per point: x, y, z, rho, ux, uy, uz)
    real x, y, z;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int ijk = k * (ny) * (nx) + j * (nx) + i;

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

    // Output in ASCII format
    std::ofstream ascii_file(ascii_filename);
    if (!ascii_file.is_open()) {
        std::cerr << "Error opening ASCII file for writing: " << ascii_filename << std::endl;
        return;
    }

    // Write header information
    ascii_file << "TITLE = \"Example: Simple 3D-Field Data\"" << std::endl;
    ascii_file << "VARIABLES = \"X\", \"Y\", \"Z\", \"RHO\", \"UX\", \"UY\", \"UZ\"" << std::endl;
    ascii_file << "ZONE I=" << nx << ", J=" << ny << ", K=" << nz << " F=POINT" << std::endl;

    // Write data
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int ijk = k * (ny) * (nx) + j * (nx) + i;

                // ascii_file << std::fixed << std::setprecision(6)  // Fixed point format with specified precision
                //            << (i + 0.5) * dx << " "
                //            << (j + 0.5) * dy << " "
                //            << (k + 0.5) * dz << " "
                //            << rhg[ijk] << " "
                //            << uxg[ijk] << " "
                //            << uyg[ijk] << " "
                //            << uzg[ijk] << std::endl;

                    ascii_file  << (i + 0.5) * dx << " "
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

void convert_all_bin_to_tec(const std::string& folder_path, real dx, real dy, real dz)
{
    for (const auto& entry : fs::directory_iterator(folder_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".bin") {
            std::string binary_filename = entry.path().string();
            std::string ascii_filename = entry.path().stem().string() + ".tec";

            // Read binary file and convert to ASCII format
            int nx, ny, nz;
            load_tecplot3D_serial_binary_to_ascii(binary_filename.c_str(), ascii_filename.c_str(), nx, ny, nz, dx, dy, dz);
        }
    }
}

int main() {
    std::string folder_path = "D:/FromTAIYI/DUGKS-PFM-GPU/100Channel/F3D/";  // Folder path containing .bin files
    real dx = 0.482188, dy = 0.482188, dz = 0.482188;  // Assume grid spacing in each direction is 1.0

    std::cout << "check " << "dx=" << dx << ", dy=" << dy << ", dz=" << dz <<std::endl

    // Call function to convert all .bin files in the folder to .tec files
    convert_all_bin_to_tec(folder_path, dx, dy, dz);

    return 0;
}
