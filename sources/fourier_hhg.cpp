#include <iostream>
#include <filesystem>
#include <chrono>

#include <nlohmann/json.hpp>

#include <mrock/utility/OutputConvenience.hpp>
#include <mrock/utility/InputFileReader.hpp>
#include <mrock/utility/better_to_string.hpp>
#include <mrock/utility/ComplexNumberIterators.hpp>

#include "FHHG/DiracSystem.hpp"
#include "FHHG/Laser/ContinuousLaser.hpp"
#include "FHHG/Laser/CosineLaser.hpp"

constexpr double target_kappa_error = 1e-3;
constexpr int n_kappa = 10;

int main(int argc, char** argv) {
    using namespace FHHG;
    using std::chrono::high_resolution_clock;

    if (argc < 2) {
		std::cerr << "Invalid number of arguments" << std::endl;//: Use mpirun -n <threads> <path_to_executable> <configfile>" << std::endl;
		return -1;
	}

    mrock::utility::InputFileReader input(argv[1]);

    /**
     * Loading configurations
     */
    const h_float temperature = input.getDouble("T");
    const h_float E_F = input.getDouble("E_F");
    const h_float v_F = input.getDouble("v_F");
    const h_float band_width = input.getDouble("band_width");
    const h_float E0 = input.getDouble("field_amplitude");
    const h_float photon_energy = input.getDouble("photon_energy");
    const std::string laser_type = input.getString("laser_type");
    const int n_z = input.getInt("n_z");
    const int max_frequency = input.getInt("max_frequency");
    const int frequency_subdivisions = input.getInt("frequency_subdivisions");
    const int n_laser_cycles = input.getInt("n_laser_cycles");

    std::unique_ptr<Laser::Laser> laser;
    if (laser_type == "continuous") {
        laser = std::make_unique<Laser::ContinuousLaser>(photon_energy, E0, max_frequency, frequency_subdivisions);
    }
    else if (laser_type == "cosine") {
        laser = std::make_unique<Laser::CosineLaser>(photon_energy, E0, max_frequency, frequency_subdivisions, n_laser_cycles);
    }
    else {
        std::cerr << "Laser type '" << laser_type << "' is not recognized!" << std::endl;
        return 1;
    }

    DiracSystem system(temperature, E_F, v_F, band_width, photon_energy);

    /**
     * Creating output dirs
     */
    auto improved_string = [](h_float number) -> std::string {
        if (std::floor(number) == number) {
            // If the number is a whole number, format it with one decimal place
            std::ostringstream out;
            out.precision(1);
            out << std::fixed << number;
            return out.str();
        }
        else {
            std::string str = mrock::utility::better_to_string(number, std::chars_format::fixed);
            // Remove trailing zeroes
            str.erase(str.find_last_not_of('0') + 1, std::string::npos);
            str.erase(str.find_last_not_of('.') + 1, std::string::npos);
            return str;
        }
    };

    const std::string BASE_DATA_DIR = "../../data/FHHG/";
    const std::string data_subdir = input.getString("data_dir")
        + "/" + laser_type + "_laser" +
        + "/T=" + improved_string(temperature)
        + "/E_F=" + improved_string(E_F)
        + "/v_F=" + improved_string(v_F)
        + "/band_width=" + improved_string(band_width)
        + "/field_amplitude=" + improved_string(E0)
        + "/photon_energy=" + improved_string(photon_energy) 
        + "/";
    const std::string output_dir = BASE_DATA_DIR + data_subdir;
    std::filesystem::create_directories(output_dir);

    /**
     * Starting calculations
     */
    high_resolution_clock::time_point begin = high_resolution_clock::now();
    std::cout << "Computing the k integrals..." << std::endl;

    std::vector<h_complex> current_density = system.compute_current_density(laser.get(), n_z, n_kappa, target_kappa_error, output_dir);
    
    high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    std::vector<h_float> current_density_real(current_density.size());
    std::vector<h_float> current_density_imag(current_density.size());

    using namespace mrock::utility;

    std::copy(make_real_part_iterator(current_density), make_real_part_iterator_end(current_density), current_density_real.begin());
    std::copy(make_imag_part_iterator(current_density), make_imag_part_iterator_end(current_density), current_density_imag.begin());

    std::vector<h_float> frequencies(laser->n_subdivisions * laser->max_frequency);
    for (int i = 0; i < frequencies.size(); ++i) {
        frequencies[i] = laser->frequency(i);
    }

    nlohmann::json data_json {
        { "time", 				                mrock::utility::time_stamp() },
        { "current_density_real",               current_density_real },
        { "current_density_imag",               current_density_imag },
        { "T",                                  temperature },
        { "E_F",                                E_F },
        { "v_F",                                v_F },
        { "band_width",                         band_width },
        { "field_amplitude",                    E0 },
        { "photon_energy",                      photon_energy },
        { "laser_type",                         laser_type },
        { "frequencies",                        frequencies },
        { "max_frequency",                      laser->max_frequency },
        { "frequency_subdivisions",             laser->n_subdivisions },
        { "n_z",                                n_z },
        { "target_kappa_error",                 target_kappa_error },
        { "n_laser_cycles",                     n_laser_cycles }
    };
    mrock::utility::saveString(data_json.dump(4), output_dir + "current_density.json.gz");

    return 0;
}