import math
import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statistics
from fractions import Fraction
import decimal
import mysql.connector as ctr
import csv

stack=list()


'''Statistics'''
def save_data(filename, data):
    with open(filename, 'w') as file:
        for item in data:
            file.write(str(item) + '\n')
def load_data(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            data.append(float(line.strip()))
    return data
def collect_ungrouped_data():
    data = []
    while True:
        try:
            value = float(input("Enter data value (or any non-numeric input to finish): "))
            data.append(value)
        except ValueError:
            break
    return data
def collect_grouped_data():
    file_path = input("Enter the path of the CSV file containing grouped data: ")
    df = pd.read_csv(file_path)
    return df
def collect_grouped_data2():
    print('x')
def calculate_mean(data):
    mean = statistics.mean(data)
    print("Mean:", mean)
    stack.append(mean)
def calculate_median(data):
    median = statistics.median(data)
    print("Median:", median)
    stack.append(median)
def calculate_mode(data):
    try:
        mode = statistics.mode(data)
        print("Mode:", mode)
    except statistics.StatisticsError as e:
        print("Mode:", str(e))
    stack.append(mode)
def calculate_std_dev(data):
    std_dev = statistics.stdev(data)
    print("Standard Deviation:", std_dev)
    stack.append(std_dev)

# Co-ordinate Geometry
def calculat1e_intercepts(equation):
    x = sp.symbols('x')
    equation = equation.replace('^', '**')  # Replace '^' with '**' for exponentiation
    expr = sp.sympify(equation)
    x_intercepts = sp.solve(expr, x)
    y_intercepts= expr.subs(x, 0)
    return x_intercepts, y_intercepts



'''Calculus'''
# Calculus Operation
def limit_at_point():
    expr = input("Enter the expression: ")
    x = sp.Symbol('x')
    point = float(input("Enter the point at which to evaluate the limit: "))
   
    limit = sp.limit(sp.sympify(expr), x, point)
    print("Limit:", limit)
    stack.append(limit)

# Calculus Operation
def check_continuity():
    expr = input("Enter the expression: ")
    x = sp.Symbol('x')
    point = float(input("Enter the point to check continuity: "))
   
    is_continuous = sp.limit(sp.sympify(expr), x, point) == sp.sympify(expr).subs(x, point)
    print("Continuity:", "Continuous" if is_continuous else "Discontinuous")

# Calculus Operation
def differentiate_expression():
    expr = input("Enter the expression: ")
    x = sp.Symbol('x')
   
    first_derivative = sp.diff(sp.sympify(expr), x)
    second_derivative = sp.diff(first_derivative, x)
    third_derivative = sp.diff(second_derivative, x)
   
    print("1st Derivative:", first_derivative)
    print("2nd Derivative:", second_derivative)
    print("3rd Derivative:", third_derivative)
    stack.append([first_derivative,second_derivative,third_derivative])

# Calculus Operation
def differentiation_at_value():
    expr = input("Enter the expression: ")
    x = sp.Symbol('x')
    value = float(input("Enter the value at which to differentiate: "))
   
    derivative = sp.diff(sp.sympify(expr), x).subs(x, value)
    print(f"Derivative at {value}:", derivative)

    stack.append(derivative)


# Calculus Operation
def indefinite_integration():
    expr = input("Enter the expression to integrate: ")
    x = sp.Symbol('x')
   
    integral = sp.integrate(sp.sympify(expr), x)
    print("Indefinite Integral:", integral)
    stack.append(integral)


# Calculus Operation
def definite_integration():
    expr = input("Enter the expression to integrate: ")
    x = sp.Symbol('x')
    a = float(input("Enter the lower limit of integration: "))
    b = float(input("Enter the upper limit of integration: "))
   
    integral = sp.integrate(sp.sympify(expr), (x, a, b))
    print("Definite Integral:", integral)
    stack.append(integral)

# Calculus Operation
def parametric_differentiation():
    t = sp.Symbol('t')
    x_expr = input("Enter x(t): ")
    y_expr = input("Enter y(t): ")
   
    x_derivative = sp.diff(sp.sympify(x_expr), t)
    y_derivative = sp.diff(sp.sympify(y_expr), t)
   
    print("x'(t):", x_derivative)
    print("y'(t):", y_derivative)
    stack.append([x_derivative,y_derivative])

# Calculus Operation
def polar_coordinates():
    r = sp.Symbol('r')
    theta_expr = input("Enter theta(r): ")
   
    r_derivative = sp.diff(sp.sympify(theta_expr), r)
   
    print("dr/dtheta:", r_derivative)
    stack.append(r_derivative)

# Calculus Operation
def find_max_min_inflexion():
    x = sp.Symbol('x')
    expr = input("Enter the expression in terms of x: ")
    first_derivative = sp.diff(sp.sympify(expr), x)
    critical_points = sp.solve(first_derivative, x)
    second_derivative = sp.diff(first_derivative, x)
    inflection_points = sp.solve(second_derivative, x)
    max_points = []
    min_points = []
    for point in critical_points:
        if second_derivative.subs(x, point) > 0:
            min_points.append(point)
        elif second_derivative.subs(x, point) < 0:
            max_points.append(point)
    print("First Derivative:", first_derivative)
    print("Second Derivative:", second_derivative)
    print("Critical Points:", critical_points)
    print("Maxima Points:", max_points)
    print("Minima Points:", min_points)
    print("Inflection Points:", inflection_points)
    stack.append([first_derivative,second_derivative,critical_points,max_points,min_points,inflection_points])

# Calculus Operation
def solve_differential_equation():
    x = sp.Symbol('x')
    y = sp.Function('y')(x)
    equation = input("Enter the differential equation (in terms of x and y): ")
    solution = sp.dsolve(sp.sympify(equation), y)
    print("Solution:", solution)
    stack.append(solution)
   

# Calculus Operation
def calculate_area_under_curve():
    x = sp.Symbol('x')
    expr = input("Enter the expression in terms of x: ")
    a = float(input("Enter the lower limit of integration: "))
    b = float(input("Enter the upper limit of integration: "))
    area = sp.integrate(sp.sympify(expr), (x, a, b))
    print("Area under Curve:", area)
    stack.append(area)

def calculus():
    print("Calculus Operations Menu:")
    print("1. Limit at a Point")
    print("2. Check Continuity")
    print("3. Differentiate Expression")
    print("4. Differentiation at Value")
    print("5. Indefinite Integration")
    print("6. Definite Integration")
    print("7. Parametric Differentiation")
    print("8. Polar Coordinates")
    print("9. Find Maxima, Minima, and Inflection Points")
    print("10. Solve Differential Equation")
    print("11. Calculate Area under Curve")
    print("0. Exit")
    choice = input("Enter your choice: ")

    if choice == '1':
        limit_at_point()
    elif choice == '2':
        check_continuity()
    elif choice == '3':
        differentiate_expression()
    elif choice == '4':
        differentiation_at_value()
    elif choice == '5':
        indefinite_integration()
    elif choice == '6':
        definite_integration()
    elif choice == '7':
        parametric_differentiation()
    elif choice == '8':
        polar_coordinates()
    elif choice == '9':
        find_max_min_inflexion()
    elif choice == '10':
        solve_differential_equation()
    elif choice == '11':
        calculate_area_under_curve()
    elif choice == '0':
        print("Exiting.")
    else:
        print("Invalid choice. Please enter a valid option.")


# Pre-calculus operation
def inverse_of_function():
    # Input the equation in terms of x and y
    equation = input("Enter the equation (in terms of x and y): ")

    # Define symbols
    x, y = sp.symbols('x y')

    # Parse the equation
    parsed_equation = sp.Eq(sp.sympify(equation), y)

    # Solve for x in terms of y (the inverse function)
    inverse_function = str(sp.solve(parsed_equation, x))

    # Print the inverse function
    print(f"Inverse Function: x = {inverse_function[1:-1]}")

# Pre-calculus operation
def law_of_sines_angle():
    angle_A_deg = float(input("Enter the known angle A in degrees: "))
    side_a = float(input("Enter the length of side a: "))
    side_b = float(input("Enter the length of side b: "))
    angle_A_rad = math.radians(angle_A_deg)
    angle_B_rad = math.asin((side_b * math.sin(angle_A_rad)) / side_a)
    angle_B_deg = math.degrees(angle_B_rad)
    print(f"The angle opposite to side b is {angle_B_deg} degrees")
    stack.append(angle_B_deg)

# Pre-calculus Operation
def law_of_sines_side():
    angle_a_deg = float(input("Enter the angle A in degrees: "))
    angle_b_deg = float(input("Enter the angle B in degrees: "))
    side_a = float(input("Enter the length of side a: "))
    angle_a_rad = math.radians(angle_a_deg)
    angle_b_rad = math.radians(angle_b_deg)
    side_b = (side_a * math.sin(angle_b_rad)) / math.sin(angle_a_rad)
    print(f"The length of side b is {side_b}")
    stack.append(side_b)

# Pre-calculus Operation
def law_of_cosines_angle():
    side_a = float(input("Enter the length of side a: "))
    side_b = float(input("Enter the length of side b: "))
    side_c = float(input("Enter the length of side c: "))
    cos_angle_c = (side_a**2 + side_b**2 - side_c**2) / (2 * side_a * side_b)
    angle_c_rad = math.acos(cos_angle_c)
    angle_c_deg = math.degrees(angle_c_rad)
    print(f"The angle opposite to side c is {angle_c_deg} degrees")
    stack.append(angle_c_deg)

# Pre-calculus Operation
def law_of_cosines_side():
    angle_c_deg = float(input("Enter the angle C in degrees: "))
    side_a = float(input("Enter the length of side a: "))
    side_b = float(input("Enter the length of side b: "))
    angle_c_rad = math.radians(angle_c_deg)
    side_c = math.sqrt(side_a**2 + side_b**2 - 2 * side_a * side_b * math.cos(angle_c_rad))
    print(f"The length of side c is {side_c}")
    stack.append(side_c)

# Pre-calculus Operation
def complex_number_operations():
    real_part_a = float(input("Enter the real part of complex number A: "))
    imaginary_part_a = float(input("Enter the imaginary part of complex number A: "))
    real_part_b = float(input("Enter the real part of complex number B: "))
    imaginary_part_b = float(input("Enter the imaginary part of complex number B: "))
    complex_a = complex(real_part_a, imaginary_part_a)
    complex_b = complex(real_part_b, imaginary_part_b)
    print("Choose an operation:")
    print("1. Addition")
    print("2. Multiplication")
    print("3. Conjugate of A")
    print("4. Conjugate of B")
    print("5. Modulus of A")
    print("6. Modulus of B")
    print('7. Powers of i')
    choice = int(input("Enter your choice: "))
    if choice == 1:
        result = complex_a + complex_b
        print(f"Addition: {result}")
        stack.append(result)
    elif choice == 2:
        result = complex_a * complex_b
        print(f"Multiplication: {result}")
        stack.append(result)
    elif choice == 3:
        result = complex_a.conjugate()
        print(f"Conjugate of A: {result}")
        stack.append(result)
    elif choice == 4:
        result = complex_b.conjugate()
        print(f"Conjugate of B: {result}")
        stack.append(result)
    elif choice == 5:
        result = abs(complex_a)
        print(f"Modulus of A: {result}")
        stack.append(result)
    elif choice == 6:
        result = abs(complex_b)
        print(f"Modulus of B: {result}")
        stack.append(result)
    elif choice == 7:
        n = int(input("Enter the power of i: "))
        result = 1j ** n
        print(f"i^{n} = {result}")
        stack.append(result)
    else:
        print("Invalid choice. Please select a valid operation.")

# Pre-calculus Operation
def vector_net_value():
    num_vectors = int(input("Enter the number of vectors: "))
    net_i = 0
    net_j = 0
    net_k = 0
    for i in range(num_vectors):
        vector_i = float(input(f"Enter the i-component of vector {i + 1}: "))
        vector_j = float(input(f"Enter the j-component of vector {i + 1}: "))
        vector_k = float(input(f"Enter the k-component of vector {i + 1}: "))
        net_i += vector_i
        net_j += vector_j
        net_k += vector_k
    print(f"The net value of the vectors is: {net_i}i + {net_j}j + {net_k}k")
    stack.append([net_i,net_j,net_k])

# Pre-calculus Operation
def vector_resultant_magnitude():
    num_vectors = int(input("Enter the number of vectors: "))
    total_x = 0
    total_y = 0
    total_z = 0
    total_magnitude = 0
    for i in range(num_vectors):
        vector_x = float(input(f"Enter the x-component of vector {i + 1}: "))
        vector_y = float(input(f"Enter the y-component of vector {i + 1}: "))
        vector_z = float(input(f"Enter the z-component of vector {i + 1}: "))
        total_x += vector_x
        total_y += vector_y
        total_z += vector_z
        total_magnitude += vector_x ** 2 + vector_y ** 2 + vector_z ** 2
    resultant_magnitude = math.sqrt(total_magnitude)
    print(f"The resultant magnitude of the vectors is: {resultant_magnitude}")
    print(f"The resultant vector is: {total_x}i + {total_y}j + {total_z}k")
    theta_x = math.degrees(math.acos(total_x / resultant_magnitude))
    theta_y = math.degrees(math.acos(total_y / resultant_magnitude))
    theta_z = math.degrees(math.acos(total_z / resultant_magnitude))
    print(f"Theta_x: {theta_x} degrees")
    print(f"Theta_y: {theta_y} degrees")
    print(f"Theta_z: {theta_z} degrees")
    stack.append([resultant_magnitude,theta_x,theta_y,theta_z])

# Pre-calculus Operation
def sequence_nth_term():
    print("Choose a sequence:")
    print("1. Arithmetic Progression (AP)")
    print("2. Geometric Progression (GP)")
    print("3. Fibonacci Sequence")
    print("4. Harmonic Progression (HP)")
   
    choice = int(input("Enter your choice: "))
   
    if choice == 1:
        a = float(input("Enter the first term (a): "))
        d = float(input("Enter the common difference (d): "))
        n = int(input("Enter the term number (n): "))
       
        nth_term = a + (n - 1) * d
        print(f"The {n}th term of the Arithmetic Progression is: {nth_term}")
        sum_of_terms = n * (2 * a + (n - 1) * d) / 2
        print(f"The sum of the first {n} terms of the Arithmetic Progression is: {sum_of_terms}")
        stack.append([nth_term,sum_of_terms])
   
    elif choice == 2:
        a = float(input("Enter the first term (a): "))
        r = float(input("Enter the common ratio (r): "))
        n = int(input("Enter the term number (n): "))
        nth_term = a * (r ** (n - 1))
        print(f"The {n}th term of the Geometric Progression is: {nth_term}")
        if r == 1:
            sum_of_terms = n * a
        else:
            sum_of_terms = a * (1 - r ** n) / (1 - r)
        print(f"The sum of the first {n} terms of the Geometric Progression is: {sum_of_terms}")
        stack.append([nth_term,sum_of_terms])
   
   
    elif choice == 3:
        n = int(input("Enter the term number (n): "))
       
        def fibonacci(n):
            if n <= 0:
                return 0
            elif n == 1:
                return 1
            else:
                return fibonacci(n - 1) + fibonacci(n - 2)
       
        nth_term = fibonacci(n)
        print(f"The {n}th term of the Fibonacci Sequence is: {nth_term}")
        sum_of_terms = sum(fibonacci(i) for i in range(1, n + 1))
        print(f"The sum of the first {n} terms of the Fibonacci Sequence is: {sum_of_terms}")
        stack.append([ nth_term ,sum_of_terms])
   
    elif choice == 4:
        a = float(input("Enter the first term (a): "))
        d = float(input("Enter the common difference (d): "))
        n = int(input("Enter the term number (n): "))
       
        nth_term = 1 / (a + (n - 1) * d)
        print(f"The {n}th term of the Harmonic Progression is: {nth_term}")
        sum_of_terms = sum(1 / (a + (i - 1) * d) for i in range(1, n + 1))
        print(f"The sum of the first {n} terms of the Harmonic Progression is: {sum_of_terms}")
        stack.append([nth_term,sum_of_terms])
   
    else:
        print("Invalid choice. Please select a valid sequence.")


# Algebra operations resolved
def multivariable_equations():
    x, y = sp.symbols('x y')
    equation1 = input("Enter equation 1 (in terms of x and y) (eg: x -2*y +5): ")
    equation2 = input("Enter equation 2 (in terms of x and y) (eg: -2*x +y -3): ")
    solution = sp.solve([sp.Eq(sp.sympify(equation1), 0), sp.Eq(sp.sympify(equation2), 0)], (x, y))
    print("Solution:", solution)  
    stack.append(solution)

# Algebra operations resolved
def linear_equations():
    x = sp.symbols('x')
    equation = input("Enter linear equation (in terms of x) (eg: 3*x -5): ")  
    solution = sp.solve(sp.Eq(sp.sympify(equation), 0), x)
    print("Solution:", solution)
    stack.append(solution)

# Algebra operations resolved
def quadratic_equations():
    x = sp.symbols('x')
    equation = input("Enter quadratic equation (in terms of x)(eg: x**2 - 4*x + 3): ")    
    solution = sp.solve(sp.Eq(sp.sympify(equation), 0), x)
    print("Solutions:", solution)
    stack.append(solution)

# Algebra operations solved
def cubic_equations():
    x = sp.symbols('x')
    equation = input("Enter cubic equation (in terms of x)(eg: x**3 - 4*x + 3): ")
    solution = sp.solve(sp.Eq(sp.sympify(equation), 0), x)
    print("Solutions:", solution)
    stack.append(solution)

# Algebra operations
def system_of_inequalities():
    x, y = sp.symbols('x y')
    inequality1 = input("Enter inequality 1 (in terms of x and y)(eg: x + y >= 3): ")
    inequality2 = input("Enter inequality 2 (in terms of x and y)(eg: 2*x - y <= 5): ")
    solution = sp.solve_univariate_inequality(sp.sympify(inequality1), x, relational=False) & sp.solve_univariate_inequality(sp.sympify(inequality2), x, relational=False)    
    print("Solution:", solution)
    stack.append(solution)

# Algebra operations
def polynomial_factorisation():
    x = sp.symbols('x')
    polynomial = input("Enter polynomial equation (in terms of x)(eg: x**2 - 4*x + 4): ")
    factors = sp.factor(sp.sympify(polynomial))
    print("Factorization:", factors)
    stack.append(factors)

def basic_operations():
    print("Select operation:")
    print("1. Addition(2)")
    print("2. Subtraction(2)")
    print("3. Multiplication(2)")
    print("4. Division(2)")
    print("5. Quotient and Remainder(2)")
    print("6. Factorial(1)")
    print("7. Permutations(2)")
    print("8. Combinations(2)")
    print("9. Greatest Integer Function(1)")
    print("10. Smallest Integer Function(1)")
    print("11. Scientific Logarithm - log(1)")
    print("12. Natural Logarithm - ln(1)")
    print("13. Exponential Function(1)")
    print("14. Convert to Scientific Notation(0)")
    print("15. Convert to Standard Notation(0)")
    print("16. HCF(2)")
    print("17. LCM(2)")
    print("18. Convert Natural Log to Scientific Log(1)")
    print("19. Convert Decimals to Fractions(0)")
    print("20. Convert Fractions to Decimals(0)")
    print("21. Square Root(1)")
    print("22. Cube Root(1)")
    print("23. Square(1)")
    print("24. Cube(1)")
    print("25. Convert Mixed to Improper Fraction(0)")
    print("26. Convert Improper to Mixed Fraction(0)")
    print("27. Number Inverse (1/x)(1)")
    print("28. Simple Interest(0)")
    print("29. Compound Interest(0)")
    choice = input("Enter choice (1-29): ")
   
    print('For single and null operators type "0" for "Enter First/Second Number:" ')
    num1 = float(input("Enter first number: "))
    num2 = float(input("Enter second number: "))

    if choice == '1':
        result = num1 + num2
        stack.append(result)
    elif choice == '2':
        result = num1 - num2
        stack.append(result)
    elif choice == '3':
        result = num1 * num2
        stack.append(result)
    elif choice == '4':
        result = num1 / num2
        stack.append(result)
    elif choice == '5':
        quotient = num1 // num2
        remainder = num1 % num2
        print("Quotient:", quotient)
        print("Remainder:", remainder)
        stack.append([quotient,remainder])
    elif choice == '6':
        result = math.factorial(int(num1))
        stack.append(result)
    elif choice == '7':
        result = math.perm(int(num1), int(num2))
        stack.append(result)
    elif choice == '8':
        result = math.comb(int(num1), int(num2))
        stack.append(result)
    elif choice == '9':
        result = math.floor(num1)
        stack.append(result)
    elif choice == '10':
        result = math.ceil(num1)
        stack.append(result)
    elif choice == '11':
        result = math.log10(num1)
        stack.append(result)
    elif choice == '12':
        result = math.log(num1)
        stack.append(result)
    elif choice == '13':
        result = math.exp(num1)
        stack.append(result)
    elif choice == '14':
# understand the function here (important)
        num3 = int(input('enter value'))
        result = "{:.3e}".format(num3)
        stack.append(result)
    elif choice == '15':
        sci = float(input("Enter number in scientific notation (e.g., 2.5e3): "))
        print(sci)
        stack.append(sci)
    elif choice == '16':
        result = math.gcd(int(num1), int(num2))
        stack.append(result)
    elif choice == '17':
# understand the function here (important)
        result = abs(int(num1) * int(num2)) // math.gcd(int(num1), int(num2))
        stack.append(result)
    elif choice == '18':
        result = num1 * 2.303
        stack.append(result)
    elif choice == '19':
        decimal = float(input("Enter decimal number: "))
# understand the function here (important)
        result = Fraction(decimal).limit_denominator()
        stack.append(result)
    elif choice == '20':
        numerator = int(input("Enter numerator: "))
        denominator = int(input("Enter denominator: "))
        result = numerator / denominator
        stack.append(result)
    elif choice == '21':
        result = math.sqrt(num1)
        stack.append(result)
    elif choice == '22':
        result = num1 ** (1/3)
        stack.append(result)
    elif choice == '23':
        result = num1 ** 2
        stack.append(result)
    elif choice == '24':
        result = num1 ** 3
        stack.append(result)
    elif choice == '25':
        whole = int(input("Enter the whole part: "))
        num1 = int(input("Enter the numerator: "))
        num2 = int(input("Enter the denominator: "))
        numerator = (whole * num2) + num1
        denominator = num2
        result = Fraction(numerator, denominator)
        stack.append(result)
    elif choice == '26':
# understand the function here (important)
        numerator = int(input("Enter the numerator: "))
        denominator = int(input("Enter the denominator: "))
        fraction = Fraction(numerator, denominator).limit_denominator()
        whole = fraction.numerator // fraction.denominator
        numerator = fraction.numerator % fraction.denominator
        result = f"{whole} '_' {numerator}/{fraction.denominator}"
        stack.append(result)
    elif choice == '27':
        result = 1 / num1
        stack.append(result)
    elif choice == '28':
        principal = float(input("Enter Principle: "))
        rate = float(input("Enter Rate: "))
        time = float(input("Enter time (in years): "))
        result = (principal * rate * time) / 100
        stack.append(result)
    elif choice == '29':
        principal = float(input("Enter Principle: "))
        rate = float(input("Enter Rate: "))
        time = float(input("Enter time (in years): "))
        n = int(input("Enter number of times interest is compounded per year: "))
        result = principal * (1 + (rate / (n * 100))) ** (n * time) - principal
        stack.append(result)
    else:
        print("Invalid choice")
    print("Result:",result)


def trigonometry_operations():
    angle = float(input("Enter angle/value (type math.pi for pi value): "))
    print("Select operation:")
    print("1. Sine (Radians)")
    print("2. Cosine (Radians)")
    print("3. Tangent (Radians)")
    print("4. Cotangent (Radians)")
    print("5. Secant (Radians)")
    print("6. Cosecant (Radians)")
    print("7. Sine (Degrees)")
    print("8. Cosine (Degrees)")
    print("9. Tangent (Degrees)")
    print("10. Cotangent (Degrees)")
    print("11. Secant (Degrees)")
    print("12. Cosecant (Degrees)")
    print("13. Degrees to Radians")
    print("14. Radians to Degrees")
    print("15. Inverse Sine")
    print("16. Inverse Cosine")
    print("17. Inverse Tangent")
    print("18. Inverse Cotangent")
    print("19. Inverse Secant")
    print("20. Inverse Cosecant")
    choice = input("Enter choice (1-20): ")

    if choice == '1':
        result = math.sin(angle)
        stack.append(result)
    elif choice == '2':
        result = math.cos(angle)
        stack.append(result)
    elif choice == '3':
        result = math.tan(angle)
        stack.append(result)
    elif choice == '4':
        result = 1 / math.tan(angle)
        stack.append(result)
    elif choice == '5':
        result = 1 / math.cos(angle)
        stack.append(result)
    elif choice == '6':
        result = 1 / math.sin(angle)
        stack.append(result)
    elif choice == '7':
        result = math.sin(math.radians(angle))
        stack.append(result)
    elif choice == '8':
        result = math.cos(math.radians(angle))
        stack.append(result)
    elif choice == '9':
        result = math.tan(math.radians(angle))
        stack.append(result)
    elif choice == '10':
        result = 1 / math.tan(math.radians(angle))
        stack.append(result)
    elif choice == '11':
        result = 1 / math.cos(math.radians(angle))
        stack.append(result)
    elif choice == '12':
        result = 1 / math.sin(math.radians(angle))
        stack.append(result)
    elif choice == '13':
        result = math.radians(angle)
        stack.append(result)
    elif choice == '14':
        result = math.degrees(angle)
        stack.append(result)
    elif choice == '15':
# understand the function here (important)
        result_degrees = math.degrees(math.asin(angle))
        result_radians = math.asin(angle)
        result = f"Degrees: {result_degrees}, Radians: {result_radians}"
        stack.append(result)
    elif choice == '16':
        result_degrees = math.degrees(math.acos(angle))
        result_radians = math.acos(angle)
        result = f"Degrees: {result_degrees}, Radians: {result_radians}"
        stack.append(result)
    elif choice == '17':
        result_degrees = math.degrees(math.atan(angle))
        result_radians = math.atan(angle)
        result = f"Degrees: {result_degrees}, Radians: {result_radians}"
        stack.append(result)
    elif choice == '18':
# understand the function here (important)
        result_degrees = math.degrees(math.atan(1 / angle))
        result_radians = math.atan(1 / angle)
        result = f"Degrees: {result_degrees}, Radians: {result_radians}"
        stack.append(result)
    elif choice == '19':
        result_degrees = math.degrees(math.acos(1 / angle))
        result_radians = math.acos(1 / angle)
        result = f"Degrees: {result_degrees}, Radians: {result_radians}"
        stack.append(result)
    elif choice == '20':
        result_degrees = math.degrees(math.asin(1 / angle))
        result_radians = math.asin(1 / angle)
        result = f"Degrees: {result_degrees}, Radians: {result_radians}"
        stack.append(result)
    else:
        print("Invalid choice")
        return
    print("Result:", result)


def pre_algebra_operations():
    print("Select operation:")
    print("1. Factors")
    print("2. Multiples")
    print("3. Ratios")
    print("4. Rates")
    print("5. Percentage")
    print("6. Direct and Inverse Relations")
    choice = input("Enter choice (1-6): ")

    if choice == '1':
        num = int(input('Enter Number: '))
        lst=list()
        for i in range(1, num + 1):
            if num % i == 0:
                lst.append(i)
        stack.append(lst)
        print(lst)
    elif choice == '2':
        num = int(input('Enter Number: '))
        lst=list()
        print(f"Multiples of {num}:")
        for i in range(1, 11):
            lst.append(num*i)
        stack.append(lst)
        print(lst)
    elif choice == '3':
        ratio1 = int(input("Enter first number: "))
        ratio2 = int(input("Enter second number: "))
        a = math.gcd(ratio1,ratio2)  
        print("Simplified Ratio:", (ratio1/a), ":", (ratio2/a))
    elif choice == '4':
        quantity = float(input("Enter quantity (in units): "))
        time = float(input("Enter time (in seconds): "))
        rate = quantity / time
        print("Rate:", rate, "units/seconds")
    elif choice == '5':
        num = float(input("Enter numerator: "))
        denom = float(input("Enter denominator: "))
        percentage = (num / denom) * 100
        print("Percentage:", percentage, "%")
    elif choice == '6':
        x = float(input("Enter x: "))
        y = float(input("Enter y: "))  
        print("Direct Proportionality (y = kx):")
        k = y / x
        print("Constant of Proportionality (k):", k)
        print("Inverse Proportionality (xy = k):")
        k_inverse = x * y
        print("Constant of Inverse Proportionality (k):", k_inverse)
    else:
        print("Invalid choice")


def algebra_operations():
    print("Select operation:")
    print("1. Multivariable Equations")
    print("2. Linear Equations")
    print("3. Quadratic Equations")
    print("4. Cubic Equations")
    print("5. System of Inequalities")
    print("6. Polynomial Factorisation")
   
    choice = input("Enter choice (1-7): ")      

    if choice == '1':
        multivariable_equations()
    elif choice == '2':
        linear_equations()
    elif choice == '3':
        quadratic_equations()
    elif choice == '4':
        cubic_equations()
    elif choice == '5':
        system_of_inequalities()
    elif choice == '6':
        polynomial_factorisation()
    else:
        print("Invalid choice")


def co_ordinate_geometry_operations():
    print("Graphing and Math Tools Menu:")
    print("1. Plot Graph of an Equation")
    print("2. Plot Point")
    print("3. Find Intersection Points of Lines")
    print("4. Calculate Slope")
    print("5. Calculate Enclosed Area")
    print("6. Calculate Distance between Two Points")
    print("7. Calculate Domain and Range")
    print("8. Calculate X and Y Intercepts")
    print("0. Exit")
    choice = int(input("Enter your choice: "))

    if choice == 0:
        print("Exiting the program.")
        return
    elif choice == 1:
        equation = input("Enter the equation (in terms of x): ")
        x = sp.symbols('x')

        x_min = float(input("Enter the minimum x value: "))
        x_max = float(input("Enter the maximum x value: "))

        x_vals = np.linspace(x_min, x_max, 400)

        y_expr = sp.sympify(equation)
        y_lambda = sp.lambdify(x, y_expr, 'numpy')

        y_vals = y_lambda(x_vals)

        plt.plot(x_vals, y_vals)
        plt.title("Graph of Equation")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.show()
   
    elif choice == 2:

        x_point = float(input("Enter x coordinate of the point: "))
        y_point = float(input("Enter y coordinate of the point: "))

        # Create a scatter plot for the point
        plt.scatter(x_point, y_point, c='red', marker='o', label='Point')

        # Plot x and y axes with origin (0,0)
        plt.axhline(0, color='black', linewidth=0.5)  # Horizontal line (x-axis)
        plt.axvline(0, color='black', linewidth=0.5)  # Vertical line (y-axis)

        # Set the axis limits to show both negative and positive sides
        plt.xlim(-max(abs(x_point), 1) - 1, max(abs(x_point), 1) + 1)
        plt.ylim(-max(abs(y_point), 1) - 1, max(abs(y_point), 1) + 1)

        # Add labels to the x and y axes
        plt.xlabel("x")
        plt.ylabel("y")

        # Add a legend
        plt.legend()

        # Display the plot
        plt.title("Plot Point with Respect to Origin")
        plt.grid(True)
        plt.show()
 
    elif choice == 3:
        # User input for line 1
        m1 = float(input("Enter slope of line 1: "))
        b1 = float(input("Enter y-intercept of line 1: "))

        # User input for line 2
        m2 = float(input("Enter slope of line 2: "))
        b2 = float(input("Enter y-intercept of line 2: "))

        # Generate x values for plotting
        x_values = np.linspace(-10, 10, 400)

        # Calculate y values for both lines
        y_values1 = m1 * x_values + b1
        y_values2 = m2 * x_values + b2

        # Set the figure size
        plt.figure(figsize=(8, 8))

        # Plot the lines
        plt.plot(x_values, y_values1, label=f'Line 1: y = {m1}x + {b1}', color='blue')
        plt.plot(x_values, y_values2, label=f'Line 2: y = {m2}x + {b2}', color='green')

        # Calculate the intersection point
        x_intercept = (b2 - b1) / (m1 - m2)
        y_intercept = m1 * x_intercept + b1

        # Plot the intersection point
        plt.plot(x_intercept, y_intercept, 'ro', label='Intersection')

        # Add labels, legend, and grid
        plt.title("Graphs of Lines and Intersection Point")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.grid(True)

        # Show the plot
        plt.show()

        # Print the coordinates of the intersection point
        print(f"Intersection Point: ({x_intercept}, {y_intercept})")

    elif choice == 4:
        x1 = float(input("Enter x coordinate of the first point: "))
        y1 = float(input("Enter y coordinate of the first point: "))
        x2 = float(input("Enter x coordinate of the second point: "))
        y2 = float(input("Enter y coordinate of the second point: "))
        slope = (y2 - y1) / (x2 - x1)
        print("Slope: ",slope)
   
    elif choice == 5:
        # User input for the number of lines
        num_lines = int(input("Enter the number of lines: "))
        slope_values = []  # List to store slope values
        y_intercept_values = []  # List to store y-intercept values

        # Collect information about each line
        for i in range(num_lines):
            m = float(input(f"Enter slope of line {i + 1}: "))
            b = float(input(f"Enter y-intercept of line {i + 1}: "))
           
            # Append slope and y-intercept values to the respective lists
            slope_values.append(m)
            y_intercept_values.append(b)

        # Calculate intersection points
        vertices = []  # List to store intersection points

        for i in range(num_lines):
            for j in range(i + 1, num_lines):
                # Extract slope and y-intercept values for the two lines
                m1 = slope_values[i]
                b1 = y_intercept_values[i]
                m2 = slope_values[j]
                b2 = y_intercept_values[j]

                # Calculate intersection point
                x_intersection = (b2 - b1) / (m1 - m2)
                y_intersection = m1 * x_intersection + b1

                # Append the intersection point to the vertices list
                vertices.append((x_intersection, y_intersection))


        # Print the vertices
        print("Vertices:", vertices)

        # Calculate the area using the Shoelace Formula
        area = 0
        n = len(vertices)
        for i in range(n):
            j = (i + 1) % n
            a = vertices[i][0]
            b = vertices[j][1]
            c = vertices[j][0]
            d = vertices[i][1]
            print(a,b,c,d)

           
            area += ((a * b) - (c * d))
            print(area)

        # Ensure the area is positive and print it
        area = 0.5 * abs(area)          
        print("Enclosed Area:", area)

    elif choice == 6:
        x1 = float(input("Enter x coordinate of the first point: "))
        y1 = float(input("Enter y coordinate of the first point: "))
        x2 = float(input("Enter x coordinate of the second point: "))
        y2 = float(input("Enter y coordinate of the second point: "))
        distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        print("Distance between the points:", distance)
   
    elif choice == 7:
        equation = input("Enter the equation (in terms of x): ")
        x_min = float(input("Enter the minimum value of x: "))
        x_max = float(input("Enter the maximum value of x: "))          
        x_vals = np.linspace(x_min, x_max, 1000)
        y_vals = [eval(equation) for x in x_vals]
        domain = [x_min, x_max]
        range_vals = [min(y_vals), max(y_vals)]           
        print("Domain:", domain)
        print("Range:", range_vals)
   
    elif choice == 8:
        # Input an equation in terms of x
        equation = input("Enter the equation (in terms of x): ")

        # Define the variable x
        x = sp.symbols('x')

        # Parse the equation
        parsed_equation = sp.sympify(equation)

        # Calculate the x-intercept (y = 0)
        x_intercept = str(sp.solve(parsed_equation, x))

        # Calculate the y-intercept (x = 0)
        y_intercept = parsed_equation.subs(x, 0)

        # Print the results
        print(f"X-intercepts: {x_intercept[1:-1]}")
        print(f"Y-intercept: {y_intercept}")
   
    else:
        print("Invalid choice. Please enter a valid option.")
        # ... (plot_equation function)

def pre_calculus_operations():
    print("1. Inverse of a function")
    print("2. Law of Sines")
    print("3. Law of Cosines")
    print("4. Complex Numbers")
    print("5. Vectors")
    print("6. Sequence and Series")
    print("0. Exit")
   
    choice = input("Enter your choice: ")
   
    if choice == '1':
        inverse_of_function()
    elif choice == '2':
        print("1. Calculate angle when one angle and two sides are given")
        print("2. Calculate side when both angles and one side are given")
        sub_choice = input("Enter sub-choice: ")
        if sub_choice == '1':
            law_of_sines_angle()
        elif sub_choice == '2':
            law_of_sines_side()
    elif choice == '3':
        print("1. Calculate angle when all three sides are given")
        print("2. Calculate opposite side when angle and two sides are given")
        sub_choice = input("Enter sub-choice: ")
        if sub_choice == '1':
            law_of_cosines_angle()
        elif sub_choice == '2':
            law_of_cosines_side()
    elif choice == '4':
        complex_number_operations()
    elif choice == '5':
        print("1. Find net value (x, y, z) of vectors")
        print("2. Calculate resultant magnitude of vectors")
        sub_choice = input("Enter sub-choice: ")
        if sub_choice == '1':
            vector_net_value()
        elif sub_choice == '2':
            vector_resultant_magnitude()
def statistics_menu():
    data = []  # This will store collected data
   
    while True:
        print("Select statistics operation:")
        print("1. Collect Ungrouped Data")
        print("2. Collect Grouped Data from CSV")
        print("3. Exit")
       
        choice = input("Enter choice (1-3): ")
       
        if choice == '1':
            data = collect_ungrouped_data()
            while True:
                print("Select operation:")
                print("1. Calculate Mean")
                print("2. Calculate Median")
                print("3. Calculate Mode")
                print("4. Calculate Standard Deviation")
                print("5. Go back to main menu")
               
                sub_choice = input("Enter sub-choice (1-5): ")
               
                if sub_choice == '1':
                    calculate_mean(data)
                elif sub_choice == '2':
                    calculate_median(data)
                elif sub_choice == '3':
                    calculate_mode(data)
                elif sub_choice == '4':
                    calculate_std_dev(data)
                elif sub_choice == '5':
                    break
                else:
                    print("Invalid choice")
        elif choice == '2':
            data = collect_grouped_data()
            while True:
                print("Select operation:")
                print("1. Calculate Mean")
                print("2. Calculate Median")
                print("3. Calculate Mode")
                print("4. Calculate Standard Deviation")
                print("5. Go back to main menu")
               
                sub_choice = input("Enter sub-choice (1-5): ")
               
                if sub_choice == '1':
                    calculate_mean(data)
                elif sub_choice == '2':
                    calculate_median(data)
                elif sub_choice == '3':
                    calculate_mode(data)
                elif sub_choice == '4':
                    calculate_std_dev(data)
                elif sub_choice == '5':
                    break
                else:
                    print("Invalid choice")
        elif choice == '3':
            return "exit"
        else:
            print("Invalid choice")


t=True
while t:
    print('******************************************')
    print("Mathematical Operations Menu:")
    print("1. Basic Operations")
    print("2. Trigonometry")
    print("3. Pre-Algebra")
    print("4. Algebra")
    print("5. Pre-Calculus")
    print("6. Co-ordinate Geometry")
    print('7. Statistics')
    print('8. Calculus')
    print("0. Exit")
    print('******************************************')

    choice = int(input("Enter your choice: "))
    if choice == 0:
        t=False

        print("Exiting the program.")
    elif choice == 1:
        basic_operations()
    elif choice == 2:
        trigonometry_operations()
    elif choice == 3:
        pre_algebra_operations()
    elif choice == 4:
        algebra_operations()
    elif choice == 5:
        pre_calculus_operations()
    elif choice == 6:
        co_ordinate_geometry_operations()
    elif choice == 7:
        statistics_menu()
    elif choice == 8:
        calculus()
    elif choice==0:
        t=False
    else:
        print('invalid')

    


