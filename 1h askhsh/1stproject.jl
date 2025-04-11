#! THIS IS THE 1st SET of Physics & Technology of Semiconductor Devices

#? PRWTI ASKISI 
# a) Eg = 1.12 eV , ni(300) = 10^10 cm^-3 for Si and - 250 Celsius =< T =< 1000 Celsius 
#where T is the temperature 

using Plots

# endogeneis sugkentrwsi ilektroniwn ni se imiagwgo dinetai apo thn sxesh :
# ni = sqrt(N_c * N_v ) * (-E_g / 2*k*T)   


#TODO HYPERPARAMETERS :
Eg = 1.12 #(eV)
ni_300 = 10^10 #(cm^-3)
k = 8.62*10^(-5) # (eV/k)
T_values = collect(23:10:1273) #with step 50 23K until 1273 Kelvin
n_i_1a =Float64[]

#! 1a) Stathero Eg 
for T in T_values
    n_i_curr = ni_300 * sqrt((T/300)^3) * exp((-Eg/(2*k))*((1/T)-(1/300)))
    push!(n_i_1a,n_i_curr)
end

#! 1b) Thermokrasiakh exartisi tou energeiakou xasmatos
#TODO HYPERPARAMETES 
Eg_0 = 1.166 # (eV)
a= 4.73*10^-4 # (eV/K)
b = 636 # (Kelvin)
n_i_1b = Float64[]

for T in T_values
    Eg_curr = Eg_0 - ((a*T^2)/(T + b))
    n_i_curr = ni_300 * sqrt(((T/300)^3)) * exp((-Eg_curr/(2*k))*((1/T)-(1/300)))
    push!(n_i_1b,n_i_curr)
end

#! 1c) Nc kai Nv analoga tou T^(3/2)
#TODO HYPERPARAMETERS
#Efoson exoume mia analigoa sta Nc kai Nv apo ton arxiko typo 
#Tha xrisimopoihsw ton idio tipo me prin apla tha exw kai to T/300 ^ 3/2
Nc_300 = 2.8 *10^19 #(cm^-3)
Nv_300 = 1*10^19 #(cm^-3)
n_i_1c = Float64[]


for T in T_values
    Eg_curr = Eg_0 - ((a*T^2)/(T + b))
    Nc = Nc_300 * (T/300)^(3/2)
    Nv = Nv_300 * (T/300)^(3/2)
    n_i_curr = sqrt(Nc * Nv) * exp(-Eg_curr/(2*k*T))
    push!(n_i_1c, n_i_curr)
end



plot_1 =plot(
    T_values, [n_i_1a, n_i_1b, n_i_1c],
    yscale = :log10,
    ylim = (1e-30,1e30),
    xlabel = "Temperature (K)",
    ylabel = "n_i (cm⁻³) *log10*",
    title = "Ενδογενής Συγκέντρωση Ηλεκτρονίων για το Si",
    labels = ["Σταθερό Eg" "Eg(T)" "Eg(T) και εξαρτηση Τ^3/2 για Nv και Nc"],
    linewidth = 3,
    color = [:darkslateblue :blue2 :limegreen],
    legend = :bottomright,
)


savefig(plot_1, "Set 1o Askhsh 1")


#? DEUTERI ASKHSH
#TODO HYPERPARAMETERS 
Nd = 10^16 #(cm^-3)
#shmatniko thetw to Ec ~= Ef afou douleyv me diafores ara Ed - Ef = Ec-0.054-Ef ara feugoun ara =-0.054
Ec = 0
Ed = Ec- 0.054 # (eV) kai to Ec einai tou pyritiou
n_ion=[] # Sygkentrwsh ionizsmenwn prosmiksewn

for T in T_values

    n_d_plus = Nd/(1+2*exp((-0.054)/(k*T)))
    push!(n_ion,n_d_plus)

end

plot_2 = plot(
    T_values,[n_ion],
    yscale = :log10,
    ylim=(1e15,1e17),
    ylabel= "Συγκέντρωση Ιονισμένων Προσμίξεων (cm^-3)",
    xlabel = "Θερμοκρασια (Τ)",
    title= "Συγκέντρωση Προσμίξεων Si με As",
    legend = :bottomleft,
    linewidth = 3,
    color = [:maroon]
)

savefig(plot_2, "Set 1o Askhsh 2.png")

#! ONE FINAL PLOT ALL OF THEM TOGETHER

plot_all1 = plot(
    T_values,[n_i_1c,n_ion],
    yscale = :log10,
    ylim=(1e0,1e20),
    ylabel= "Συγκεντρωση Ηλεκτρονιων (cm^-3)",
    xlabel = "Θερμοκρασια (Τ)",
    title = "Συγκεντρωση Ηλεκτρονιων Ενδογενων/Προσμιξεων",
    legend = :bottomleft,
    linewidth = 3,
    color = [:limegreen :maroon],
    labels = ["Ενδογενη Συγκεντρωση απο 1γ)" "Συγκεντρωση Προσμίξεων"],
)

savefig(plot_all1, "Set 1o Askhsh 2 με κανονικη κλιμακα.png")

#! ΚΑΙ ΑΛΛΙ ΜΙΑ ΜΕ 1/Τ 

plot_all2 = plot(
    1 ./T_values,[n_i_1c,n_ion],
    yscale = :log10,
    ylim = (1e-30,1e20),
    ylabel= "Συγκεντρωση Ηλεκτρονιων (cm^-3)",
    xlabel = "Θερμοκρασια (1/Τ)",
    title = "Συγκεντρωση Ηλεκτρονιων Ενδογενων/Προσμιξεων",
    legend = :bottomleft,
    linewidth = 3,
    color = [:green3 :brown],
    labels = ["Ενδογενη Συγκεντρωση απο 1γ)" "Συγκεντρωση Προσμίξεων"],
)

savefig(plot_all2,"Set 1o Askhsh 2 με T^-1 κλιμακα.png")
