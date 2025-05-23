"""
This is the second set of excercises for the course Φυσικη και Τεχνολογια Ημιαγωγικων Διαταξεων

"""
#! 2h ERGASIA
# This julia file is where the analysis and plotting will happend :

using Plots
using DataFrames
using XLSX

#! ------------- 2h Erwtisi ------------------
#? a erwtima --------->
#diatomi
A = 0.04 * 0.04 #cm * cm

# Main Dataframe 

# Anastrofi polwsi V_r 
V_r = [0.0, 0.5, 1.0, 1.5, 2.0 ,3.0, 4.0, 5.0]

# Xwritikothta
C = [3.52e-11, 2.68e-11, 2.24e-11, 1.98e-11, 1.81e-11, 1.55e-11, 1.35e-11, 1.25e-11]

# Διηλεκτρικη σταθερα 
k = 12

# Bohthitikos syntelestis
boi_syt = Float64[]

for i in eachindex(C) 
    push!(boi_syt, 1/(C[i])^2)
end
    
print(boi_syt)

# Ypologismos a :

sum_x = sum(V_r)
sum_y = sum(boi_syt)
sum_yx = sum((boi_syt) .* V_r)
sum_x2 = sum(V_r .* V_r)
n = length(V_r)

a = ((n * sum_yx) - (sum_x * sum_y))/(n * sum_x2 - (sum_x)^2) 

# Ypologismos b :
b = (sum_y - a * sum_x)/n 
#thimizw a kai b exoun syntelests F^-2*V^-1 kai F^-2 antistoixws

# Ypologismos Vbi :
Vbi = - b/a 

# Ypologismos a kai b apo grammikh palindromisi
# y = 1/C^2 kai x = V_r

#? PLOTTING
fit_line = @. a * V_r + b

plot_2a = scatter(V_r, boi_syt,
    xlabel = "Ανάστροφη Πόλωση V_r (V)",
    ylabel = "1/C² (F⁻²)",
    label = "Δεδομένα",
    title = "Διάγραμμα 1/C² - V_r με παλινδρόμηση",
    marker = :hexagon,
    color = :orange2)

plot!(V_r, fit_line,
    label = "Γραμμή παλινδρόμησης",
    linewidth = 3,
    color = :blueviolet)


# Αποθήκευση
savefig(plot_2a, "../Fysikh-kai-Texnologia-Imiagwgwn/2h askhsh/plots/2Α_Ερωτημα.png")

#? β Ερωτημα ---------------->
# Να = Νd = Ν efosonexoume idio arithmo prosmijewn n kai p

q = 1.6 * 10e-19 #Fortio metrimeno se [Coulomb]
e_o = 8.854e-14 #Dihlektriki stathera sto keno [F/cm]
e_s = k * e_o #F/cm

N = 2/(a * q * e_s * A^2) #Mia taji megethous katw apo thn endeiktikh apanthsh 10^14  kai stis duo perioxes


#? γ Ερωτημα ---------------->
#Kata metro tha exoyme idia Q- kai Q+

Q = Float64[]

V_r
Vbi
C
C[2]

for i in eachindex(C)
    push!(Q, C[i] *(Vbi + V_r[i]) ) 
end

Q

print(Q)

#? δ Ερωτημα ------------------->

# prwta ypologsimos toy W 
W = Float64[]

for i in eachindex(C)
    push!(W, (e_s*A)/C[i])
end 

W

#kai epeita to Emax 

E_max = Float64[]

for i in eachindex(C)
    push!(E_max, (Vbi + V_r[i])/W[i])
end

E_max

#? ε Ερωτημα -------------------->

V_x_0 = Float64[]

for i in eachindex(C)
    push!(V_x_0, -(q*N*W[i])/(2*e_s))
end

V_x_0

#? στ Ερωτημα ---------------------------> ΠΛΟΤΣ
   
#γ Ερωτημα ->

plot_2c = plot(
    V_r, Q,
    xlabel = "Ανάστροφη Πολωση V_r (V)",
    ylabel = "|Q|",
    title = "Φορτιο για καθε v_r",
    linewidth = 3,
    marker = :hexagon,
    markercolor = :orange2,
    color = :blueviolet
)   

savefig(plot_2c, "../Fysikh-kai-Texnologia-Imiagwgwn/2h askhsh/plots/2Γ_Ερωτημα.png")

# δ Ερωτημα ->

plot_2d_1 = plot(
    V_r, W,
    xlabel = "Ανάστροφη Πολωση V_r (V)",
    ylabel = "Πλατος περιοχης Απογυμνωσης",
    title = " Πλατος περιοχης Απογυμνωσης για καθε v_r",
    linewidth = 3,
    marker = :hexagon,
    markercolor = :orange2,
    color = :blueviolet
)

savefig(plot_2d_1, "../Fysikh-kai-Texnologia-Imiagwgwn/2h askhsh/plots/2Δ_Ερωτημα_1.png")

plot_2d_2 = plot(
    V_r, E_max,
    xlabel = "Ανάστροφη Πολωση V_r (V)",
    ylabel = "Μεγιστο Ηλεκτρικο πεδιο",
    title = " Μεγιστο Ηλεκτρικο Πεδιο για καθε v_r",
    linewidth = 3,
    marker = :hexagon,
    markercolor = :orange2,
    color = :blueviolet
)

savefig(plot_2d_2, "../Fysikh-kai-Texnologia-Imiagwgwn/2h askhsh/plots/2Δ Ερωτημα 2.png")

# ε Ερωτημα ->

plot_2e = plot(
    V_r, V_x_0,
    xlabel = "Ανάστροφη Πολωση V_r (V)",
    ylabel = "Δυναμικο στην μεταλλουργικη Επαφη",
    title = " Δυναμικο στην μεταλλουργικη Επαφη για καθε v_r",
    linewidth = 3,
    marker = :hexagon,
    markercolor = :orange2,
    color = :blueviolet
)

savefig(plot_2e, "../Fysikh-kai-Texnologia-Imiagwgwn/2h askhsh/plots/2E Ερωτημα.png")

#!-------------------------------- 3h Erwtisi -----------------------------------------

#? α) Ερωτημα ------------------------------------>

C_0 = 8e-12 #F
C_5 = 3.2e-12 #F 

#boithitioi syntelestes

boi_syt_0 = 1/(C_0^2)
boi_syt_5 = 1/(C_5^2)

#Euresi klishs a kai tomis b

a_ = (boi_syt_5 - boi_syt_0)/(-5-0)

b_ = (boi_syt_0)

V_bi = -b_/a_ #einai konta stin timi 1 eV

#? β) Ερωτημα -------------------------------------->
A = 1.6e-4 #cm^2
N_L = 2/(q*e_s*A^2*abs(a_)) #cm^-3

#? γ) Ερωτημα ----------------------------------->
T = 300 # K 
k =  1.38e-23 # 
q = 1.6e-19
n_i = 1.5e10 #cm^-3
V_bi
idk = k*T/q

N_H = ((n_i^2) / (N_L))*exp(V_bi/idk) #3 tajeis megethous panw 


#! --------------------------------- 1h Erwtisi (BONUS) ---------------------------------------

xf = XLSX.openxlsx("./Fysikh-kai-Texnologia-Imiagwgwn/2h askhsh/double-diode2025-start example.xlsx")

sh = xf["ddiode"]
sh[:]

data = XLSX.getdata(sh)
data
Voltage = data[2:46] #Voltages (Volt)
Currents = data[2:46,2] # Currents (A)

ln_I = Float64[]

for i in eachindex(Currents)
    push!(ln_I,log(Currents[i]))
end

println(ln_I)

plot_1_1 = plot(
    Voltage,ln_I,
    xlabel = "Voltage of the Diode (V)",
    ylabel = "Logarithmic Current of the diode ln(I)",
    title = " ln(I) - Voltage",
    linewidth = 5,
    marker = :hexagon,
    markercolor = :orange2,
    color = :blueviolet
)

savefig(plot_1_1, "./Fysikh-kai-Texnologia-Imiagwgwn/2h askhsh/plots/1 ln(I) V.png")

I_s_low = Float64[]
I_s_high = Float64[]
n_l = 1
n_h = 2

for i in eachindex(Voltage)
    push!(I_s_low, Currents[i]/(exp(Voltage[i]/(n_l*idk))-1))
    push!(I_s_high, Currents[i]/exp(Voltage[i]/(n_h*idk)))
end

I_s_low
I_s_high