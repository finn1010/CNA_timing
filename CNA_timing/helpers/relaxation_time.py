from sympy import symbols, Matrix, Identity, det, solve

mu, gamma, lambd = symbols('mu gamma lambda')

RateMatrix = Matrix([
    [-4*mu, gamma, 0, 0, 0],                
    [4*mu, -(3 * mu + gamma), 2 * gamma, 0, 0], 
    [0, 3 * mu, -(2 * mu + 2*gamma), 3*gamma, 0], 
    [0, 0, 2*mu, -(3*gamma+mu), 4*gamma],
    [0, 0, 0, mu, -4*gamma]           
])

char_eq = det(RateMatrix - lambd * Identity(5))

eigenvalues = solve(char_eq, lambd)

eigenvalues_symbolic = RateMatrix.eigenvals()

print("Characteristic Equation:")
print(char_eq)
print("\nEigenvalues (explicit solution):")
print(eigenvalues)
print("\nEigenvalues (symbolic computation):")
print(eigenvalues_symbolic)
