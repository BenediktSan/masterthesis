all: 	
	mpicxx main.cpp  -o main && mpirun -n 4 main
	python ../python_files/single_system_plots.py
	python ../python_files/single_system_zeeman_plots.py
	python ../python_files/averaged_system_correlation_plots.py
	python ../python_files/averaged_system_pair_correlation_plots.py
	python ../python_files/comparison_plot.py
	python ../python_files/correlation_convergence.py




build/plots/plot1.pdf:	build/daten/Magnetization.txt	\ |build\
	python test_plots.py

build/daten/Magnetization.txt: \
	g++ main.cpp  -o main && "home/sander/Data/Masterarbeit/main/"main


run_all:
	g++ main.cpp  -o main && "/home/sander/Data/Masterarbeit/main/"main
	python ../python_files/single_system_plots.py
	python ../python_files/single_system_zeeman_plots.py
	python ../python_files/averaged_system_plots.py

run_testfile:
	g++ testfile.cpp -o testfile && "/home/sander/Data/Masterarbeit/main/"testfile

run_mpi:
	mpicxx main.cpp  -o main && mpirun -n 4 main "3D"
	python ../python_files/single_system_plots.py > /dev/null
	python ../python_files/single_system_zeeman_plots.py > /dev/null
	python ../python_files/averaged_system_correlation_plots.py > /dev/null 
	python ../python_files/averaged_system_pair_correlation_plots.py > /dev/null
	python ../python_files/comparison_plot.py > /dev/null
	python ../python_files/correlation_convergence.py > /dev/null
	python ../python_files/variance_for_certain_step.py > /dev/null
	python ../python_files/pair_correlation_comparison.py > /dev/null
	python ../python_files/compare_certain_pair_correlations.py > /dev/null
	python ../python_files/pair_correlation_convergence.py > /dev/null
	python ../python_files/hybrid_plots.py > /dev/null


run_mpi_2D:
	mpicxx main.cpp  -o main && mpirun -n 4 main "2D"
	python ../python_files/2D_plots/averaged_system_correlation_plots_2D.py > /dev/null
	python ../python_files/2D_plots/averaged_system_pair_correlation_plots_2D.py > /dev/null
	python ../python_files/2D_plots/pair_correlation_convergence_2D.py > /dev/null
	python ../python_files/2D_plots/variance_for_certain_step_2D.py > /dev/null
	python ../python_files/2D_plots/pair_correlation_comparison_2D.py > /dev/null
	python ../python_files/2D_plots/compare_certain_pair_correlations_2D.py > /dev/null
	python ../python_files/2D_plots/pair_correlation_comparison_2D_external_data.py > /dev/null
	python ../python_files/2D_plots/pair_correlation_comparison_2D_XXZ_external_data.py > /dev/null 

run_2D_single_core:
	mpicxx main.cpp  -o main && mpirun -n 1 main "2D"

run_3D_single_core:
	mpicxx main.cpp  -o main && mpirun -n 1 main "3D"


plots_3D:
	python ../python_files/single_system_plots.py
	python ../python_files/single_system_zeeman_plots.py
	python ../python_files/averaged_system_correlation_plots.py
	python ../python_files/averaged_system_pair_correlation_plots.py
	python ../python_files/comparison_plot.py
	python ../python_files/correlation_convergence.py
	python ../python_files/variance_for_certain_step.py
	python ../python_files/pair_correlation_comparison.py
	python ../python_files/compare_certain_pair_correlations.py
	python ../python_files/pair_correlation_convergence.py
	python ../python_files/hybrid_plots.py

plots_2D:
	python ../python_files/2D_plots/averaged_system_correlation_plots_2D.py
	python ../python_files/2D_plots/averaged_system_pair_correlation_plots_2D.py
	python ../python_files/2D_plots/pair_correlation_convergence_2D.py
	python ../python_files/2D_plots/variance_for_certain_step_2D.py
	python ../python_files/2D_plots/pair_correlation_comparison_2D.py
	python ../python_files/2D_plots/compare_certain_pair_correlations_2D.py
	python ../python_files/2D_plots/pair_correlation_comparison_2D_external_data.py
	python ../python_files/2D_plots/pair_correlation_comparison_2D_XXZ_external_data.py

plots: plots_2D plots_3D

renew_all_plots_3D: 

	python ../python_files/averaged_system_correlation_plots.py renew_all_plots
	python ../python_files/averaged_system_pair_correlation_plots.py renew_all_plots
	python ../python_files/comparison_plot.py renew_all_plots
	python ../python_files/correlation_convergence.py renew_all_plots
	python ../python_files/variance_for_certain_step.py renew_all_plots
	python ../python_files/pair_correlation_comparison.py renew_all_plots
	python ../python_files/compare_certain_pair_correlations.py renew_all_plots
	python ../python_files/pair_correlation_convergence.py renew_all_plots
	python ../python_files/hybrid_plots.py renew_all_plots

renew_all_plots_2D:
	python ../python_files/2D_plots/averaged_system_correlation_plots_2D.py renew_all_plots
	python ../python_files/2D_plots/averaged_system_pair_correlation_plots_2D.py renew_all_plots
	python ../python_files/2D_plots/pair_correlation_convergence_2D.py renew_all_plots
	python ../python_files/2D_plots/variance_for_certain_step_2D.py renew_all_plots
	python ../python_files/2D_plots/pair_correlation_comparison_2D.py renew_all_plots
	python ../python_files/2D_plots/compare_certain_pair_correlations_2D.py renew_all_plots
	python ../python_files/2D_plots/pair_correlation_comparison_2D_external_data.py renew_all_plots
	python ../python_files/2D_plots/pair_correlation_comparison_2D_XXZ_external_data.py renew_all_plots
	

renew_ALL: renew_all_plots_3D renew_all_plots_2D

copy_from_thesis:
	python ../python_files/structural_python_files/copy_txt_files.py ../masterthesis/images/plots build/daten build/daten_2D
	python ../python_files/structural_python_files/copy_txt_files.py ../masterthesis/images/plots_2D build/daten build/daten_2D

delete_copied_files:
	python ../python_files/structural_python_files/delete_copied_files.py

renew_in_thesis: copy_from_thesis renew_ALL delete_copied_files
	printf "\e[1;32m \n\n\t\tALL PLOTS FOR MASTERTHESIS HAVE BEEN RENEWED\e[0m\n"

build:
	mkdir -p build
	mkdir -p build/daten
	mkdir -p build/daten_2D


clean:
	rm -rf build
	rm -rf copied_files.json

del_plots:
	rm -rf build/plots
	rm -rf build/plots2D

FORCE:

.PHONY: all clean
