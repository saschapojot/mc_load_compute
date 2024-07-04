For running computation manually:
1. Set the potential function, create the input parameter file, and make sure
that the potential function can read the parameters.
2. set the potential function by using
   auto mcObj = mc_computation(T, std::make_shared<funcName>("rowName"), observableName);
where funcName is the potential function's name, and rowName is the row's name in the parameter file.
3. make run_mc_load_and_compute
4. ./run_mc_load_and_compute T U 
5. If you want to check observables for all T, cd checkObservables/, 
for example to check L, use: python checkL.py funcName rowName
6. If you want to check observables for one T, cd oneTCheckObservables/,
 for example to check y0, use: python checky0OneT.py funcName rowName T
7. To extract the effective data from U, L, y0,z0,y1, first: cd data2Json/
then: python dist_data2json.py funcName rowName
8. To make a plot for all temperatures, first: cd plt/
then python json2plt.py funcName rowName



For running automatically:
1. Set the potential function, create the input parameter file, and make sure
   that the potential function can read the parameters.
2. set the potential function by using
   auto mcObj = mc_computation(T, std::make_shared<funcName>("rowName"), observableName);
3. make run_mc_load_and_compute
4. python launch_mc_oneT.py T
5. To extract the effective data from U, L, y0,z0,y1, first: cd data2Json/
   then: python dist_data2json.py funcName rowName
6. To make a plot for all temperatures, first: cd plt/
   then python json2plt.py funcName rowName