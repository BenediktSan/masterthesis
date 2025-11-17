#size of labels and ticks in plots
size_label = 17

#colors for plots
exp_color = "black"
FID_color = "royalblue"
timo_color_FID =  "forestgreen"
starkov_color_FID = "red"
time_exp_color = "orange"
przemek_color = "purple"

fit_color = "gray"

#hybrid colors
optimal_FID = "darkorange"
suboptimal_FID_1 = "red"
suboptimal_FID_2 = "green"
hybrid_deviation = FID_color

######pair corr colors

pair_corr_sum_color = "gold"
FID_color_pair_corr = "black"

pair_color_2D = FID_color



#labels for the axis
time_label = r'$t$ / $µs$'
FID_label = r'$\mathcal{F}(t)$ '


#compare.py
inverse_system_size_label = r'$1/L^3$ '
deviation_3D_label = r'$\Delta \mathcal{F}(t)$ / $\%$ '


#convergence.py
interim_amount_label = r'$N$'
interim_amount_label_inverse = r'$1/N$'
deviation_3D_label_interim = deviation_3D_label

#hybrid.py
deviation_3D_label_hybrid = r'Deviation in $\%$ '
truncation_label = r'Cut-off index $c_\mathrm{trunc}$'


#pair corr comparison labels    

pair_correlation_label = r'$m_c$ $G_c^{xx}(t)$ '
certain_pair_correlation_label = r'$G_c^{xx}(t)$ '

deviation_3D_label_pair_corr = r'$\Delta G_c^{xx}(t)$ / $\%$ '





######## 2D files

time_label_2D = r'$t$ / $J^{-1}$'
pair_corr_label_2D = pair_correlation_label
FID_label_2D = FID_label

deviation_2D_label_interim = deviation_3D_label_interim



#magentic field labels
B001 = r'$ {B}= (0,0,1)^T $'
B011 = r'$ {B}= (0,1,1)^T $'
B111 = r'$ {B} = (1,1,1)^T $'




#some functions for plots

#function to gain ticks for inverse system size plots

def gain_ticks_inverse_system_size(particle_amount):
    ticks = []
    tick_labels = []
    for i in range(0, len(particle_amount)):
        ticks.append( 1/(particle_amount[i]**3) )
        number_string = str(particle_amount[i]) + r"$^3$"
        label_string = r"$1/$" + number_string
        tick_labels.append( label_string )
    return ticks, tick_labels