#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <omp.h>
#include <chrono>

#include "Density.h"
#include "DensityG.h"
#include "DFT.h"
#include "pk.h"

int main(int argc, char* argv[])
{
    //const char* argv[6] = { "valami", "configuration3d.dat", "8", "50", "64", "c"};
    //int argc = sizeof(argv) / sizeof(argv[0]);
    if (argc < 6)
    {
        std::cout << "Not enough parameters!" << std::endl;
        exit(-1);
    }
    else
    {
        std::string filename = argv[1];
        int numberOfThreads = atoi(argv[2]);
        double side = atoi(argv[3]);
        int sidelength = atoi(argv[4]);
        std::string type = argv[5];

        if (type != "c" && type != "g")
        {
            std::cout << "Type is invalid!" << std::endl;
            exit(-1);
        }

        std::ifstream file(filename);
        std::string line = "";

        int numberOfLines = 1;
        std::getline(file, line);
        int dim = std::count(line.begin(), line.end(), ' ') + 1;
        while (std::getline(file,line))
            ++numberOfLines;


        file.close();
        file.open(filename);

        DFT dft(dim, sidelength);

        std::vector<std::vector<double>> c(numberOfLines, std::vector<double>(dim));
        for (size_t i = 0; i < numberOfLines; ++i)
        {
            for (size_t j = 0; j < dim; ++j)
            {
                file >> c[i][j];
            }
        }
        file.close();

        //Density d(dim, sidelength, side);
        //d.EvaluateDensity(c.begin(), c.end());

        //std::ofstream save("results"+std::to_string(dim) + "d.dat");

        //for (double x : d.GetDensity())
        //{
        //    save << x << std::endl;
        //}
        //save.close();

        std::vector<double> density = Density(dim, sidelength, side).GetDensity();

        /*auto start = std::chrono::high_resolution_clock::now();*/

        if (type == "c")
        {
            #pragma omp parallel num_threads(numberOfThreads)
            {
                int tid = omp_get_thread_num();
                int nthreads = omp_get_num_threads();

                int chunk = c.size() / nthreads;
                int start = tid * chunk;
                int end = (tid == nthreads - 1) ? c.size() : start + chunk;

                Density d(dim, sidelength, side);
                d.EvaluateDensity(c.begin() + start, c.begin() + end);
                #pragma omp critical
                {
                    std::transform(density.begin(), density.end(), d.GetDensity().begin(), density.begin(), [](double d1, double d2) {return d1 + d2; });
                }
            }
        }
        else
        {
            #pragma omp parallel num_threads(numberOfThreads)
            {
                int tid = omp_get_thread_num();
                int nthreads = omp_get_num_threads();

                int chunk = c.size() / nthreads;
                int start = tid * chunk;
                int end = (tid == nthreads - 1) ? c.size() : start + chunk;

                DensityG d(dim, sidelength, side, 1e-7);
                d.EvaluateDensity(c.begin() + start, c.begin() + end);
                #pragma omp critical
                {
                    std::transform(density.begin(), density.end(), d.GetDensity().begin(), density.begin(), [](double d1, double d2) {return d1 + d2; });
                }
            }
        }

        double averageDensityInv = density.size() / std::accumulate(density.begin(), density.end(), 0.);
        std::transform(density.begin(), density.end(), density.begin(), [&](double d) {return d * averageDensityInv - 1.; });

        std::vector<comp> dft_input(density.size());
        std::transform(density.begin(), density.end(), dft_input.begin(), [](double d) {return comp(d, 0.);});

        std::vector<comp> dft_output = dft.FourierTransform(dft_input);

        std::string savename = "powspec" + std::to_string(dim) + "d" + std::to_string(numberOfThreads); // formerly here I wrote "results" instead of "powspec", when I tested the particle mesh at the 3. biweek. At the bottom of this file, trace of saving dft_output from the 4. biweek can also be found.
        if (type == "g")
        {
            savename += "G";
        }
        else
        {
            savename += "C";
        }

        //std::ofstream save(savename + ".dat");

        //for (comp x : dft_output)
        //{
        //    save << std::setprecision(16) << x.real() << " " << x.imag() << std::endl;
        //}

        std::vector<double> pk_grid(dft_output.size());
        double factor = std::pow(side / (sidelength * sidelength), 3);
        std::transform(dft_output.begin(), dft_output.end(), pk_grid.begin(), [&](comp d) {return factor * (d.real() * d.real() + d.imag() * d.imag());});
        std::vector<double> kmod(fourier_grid_box_3d(sidelength, side, dim));
        std::vector<double> k_centres, pk_binned;
        std::tie(k_centres, pk_binned) = bin_isotropic_modes(pk_grid, kmod, sidelength, side);
        SaveVectors(savename + ".txt", k_centres, pk_binned);

        //auto stop = std::chrono::high_resolution_clock::now();

        //std::cout << std::endl;
        //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << std::endl;

        //std::string savename = "dft_results" + std::to_string(dim) + "d" + std::to_string(numberOfThreads);
        //if (type == "g")
        //{
        //    savename += "G";
        //}
        //else
        //{
        //    savename += "C";
        //}

        //std::ofstream save(savename + ".dat");

        //for (comp x : dft_output)
        //{
        //    save << x.real() << " " << x.imag() << std::endl;
        //}
        //save.close();


    }

    return 0;
}