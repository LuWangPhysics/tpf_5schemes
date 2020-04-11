%-------------------------------------------------------------------------
                      MIT License
%-------------------------------------------------------------------------

Copyright (c) [2019] [Lu Wang] 
Email(no space): lu(dot)wangphysics(at)gmail(dot)com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


%-------------------------------------------------------------------------
                   Introduction of using the code
%-------------------------------------------------------------------------

This code serves as the numerical package for the paper:
"Tilted-pulse-front schemes for terahertz generation" 
https://arxiv.org/abs/2001.09291

where 5 different tilted-pulse-front schemes for terahertz generation problem 
are calculated. For a calculation job, the calculation time varys from 4 hours 
to 4 days (DESY max-well cluster, CPUs-32 Memory-256GB Model-E5-2640 v3 @ 2.60GHz 
Constraints-INTEL, V3, E5-2640, 256G) depending on the parameter settings. 

(1)It is valid for nonlinear material Lithium Niobate (LN).
(2)It is only valid for room temperature 300k. 
However, by adding the material properties in:
"/my_material" and "/my_discritization/f_t_discrit2D.m"
other nonlinear materials and temperature dependence of the nonlinear materials
can be easily implemented.
(3)The default setting is at 1018nm pump pulse. The default dispersion function 
of the pump should be valid from 500nm-3000nm or even more. However, if other
dispersion is needed, one can change the function 
"/my_material/Li_Nb_o3/cLN_ir.m"
(4)The parameters of the pump pulse can be modified in the "tpf2D_public.m" file.
(5)The parameters of the configuration of each scheme, e.g imaging system, 
 echelon step sizes, can be modified in "/my_tpf_setup/tpf_..." file correspondingly. 
The default imaging system is chosen such that the image of the grating (echelon)
overlaps with the pulse-front-tilt. 


%-------------------------------------------------------------------------
                   Parameters definition
%-------------------------------------------------------------------------

The input electric field is defined as: 
tau=tau_fwhm/sqrt(2*log(2));
E=E_0*exp(-t^2/tau^2)exp(-x^2/(2*sigma_in^2))exp(-y^2/(2*sigma_homo^2));

---------------------------------
Change the calculation length:
---------------------------------
In the my_input class ("my_input.m"), the variable "obj.d_beam_peak" controls the calculation length.
---------------------------------
Picking up the old calculation:
---------------------------------
Checkpoint is included. If a calculation is aborted for some reason, it will recognize the old checkpoint and pick up from where the checkpoint is saved automatically when launched next time.  
-If all the variables used in defining the saving name (file_name at line 55 in the "tpf2D_public.m" ) are the same, the program will recognize them as the same process and will pick up from the old checkpoint.
-however, if the entire run is finished (not aborted in between), the checkpoint will be deleted automatically. 
-If you wish to keep the checkpoint after the program finishes, uncomment the line 111 in the "tpf2D_public.m" function. 


%-------------------------------------------------------------------------
                   Check the validation of the parameter setting
%-------------------------------------------------------------------------

Since different users might be interested in different parameter settings, 
it is suggested to check the following aspects to avoid numerical erros.

(1)Check the output files
"/my_output/...._ir" and "/my_output/...._ir_fluence".
These two files plot the pump spectrum and the fluence respectively. 
Make sure that the numerical domains are large enough to cover the 
domain of the users' interest. If not, the frequency and the spatial domains
can be modified at "/my_discritization/x_discritization.m" and "/my_struct/input_Lu.m"
respectively.
 




