from Toolbox import MP
from Toolbox.Laminate import LaminateBuilder



laminateshear = LaminateBuilder([90,-45, 90, 45, 90, 0], True, True, 1)
# laminateshear = LaminateBuilder([45], False, False, 1)
print(laminateshear.D_matrix)
D_matrix = laminateshear.D_matrix
D11 = D_matrix[0, 0]
D12 = D_matrix[0, 1]
D22 = D_matrix[1, 1]
D66 = D_matrix[2, 2]

a = 450
b = 150
Nxycrit = 27.396 * (b/a**3) * (D11 + 2*(D12 + 2*D66) * ((a^2)/(b^2)) + D22*((a^4)/(b^4)))
print('Nxycrit = ', Nxycrit)

weight = laminateshear.h * a * b * 3 * MP.rho
print('weight of plate: ', weight, ' grams')
laminateshear.Loads = [0, 0, -10, 0 ,0 ,0]
if any(laminateshear.FailureAnalysis()[0]) >= 1:
    print('first ply failure')

EI_min = 225 * (
        (1 / 3) * (D11 + 18 * (D12 + 2 * D66) + 81 * D22)
        - (D11 + 2 * (D12 + 2 * D66) + D22)
)

print(EI_min)
print(EI_min/142000)

laminateflange = LaminateBuilder([45,-45, 0, 0, 0], True, True, 1)
D66 = laminateflange.D_matrix[2, 2]

t = 1.35 #mm
b = 30
Ftot = 10000 # newton
Nx = Ftot/b
# Analyse FPF for flange:
laminateflange.Loads = [1.9*Nx, 0, 0, 0, 0, 0]
print(laminateflange.FailureAnalysis()[0])

N_xcrit = 12 * ((b / t) ** 0.717 * D66) / (1.63 * b ** 2)
Ftot = N_xcrit*b
print('crippling final load:', Ftot)
laminateflange.CalculateEquivalentProperties()
Ex = laminateflange.Ex

EI = Ex*(1.35*30**3)/12
print(EI)

