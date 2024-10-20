import numpy as np
from sympy import symbols, lambdify

# Fungsi untuk perhitungan Gauss Jordan
def gaussJordan():
    while True:  # Perulangan untuk memungkinkan perhitungan berulang
        n = int(input("Masukkan jumlah baris dan kolom matriks A: "))

        A = []
        for i in range(n):
            row = []
            for j in range(n + 1):
                row.append(float(input(f"Masukkan elemen A[{i + 1}][{j + 1}]: ")))
            A.append(row)

        A = np.array(A)

        print("Matriks A")
        print("A = ", A)

        for i in range(n):
            max_row = i
            for k in range(i + 1, n):
                if abs(A[k][i]) > abs(A[max_row][i]):
                    max_row = k
            A[[i, max_row]] = A[[max_row, i]]

            pivot = A[i][i]
            A[i] = A[i] / pivot

            for j in range(n):
                if j != i:
                    factor = A[j][i]
                    A[j] = A[j] - factor * A[i]

        print("\nMatriks Echelon Form (Gauss-Jordan)")
        print(A)

        print("\nSolusi:")
        for i in range(n):
            print(f"x[{i + 1}] = {A[i][-1]}")

        choice = input("Lakukan perhitungan Gauss-Jordan lagi? (y/n): ")
        if choice.lower() != "y":
            break  # Keluar dari perulangan jika pengguna tidak ingin melanjutkan

# Fungsi untuk perhitungan Biseksi
# Fungsi untuk perhitungan Biseksi
def biseksi():
    def f(x, equation):
        try:
            result = eval(equation)
            return result
        except Exception as e:
            print(f"Terjadi kesalahan dalam evaluasi fungsi: {e}")
            return None  # Mengembalikan None jika terjadi kesalahan dalam evaluasi persamaan.

    while True:  # Perulangan tambahan untuk menu Biseksi
        print('Masukkan persamaan f(x). Misal: x**3-x**2-4 ')
        equation = input('(x sebagai variabel): ')
        x0 = float(input('Tebakan pertama    : '))
        x1 = float(input('Tebakan kedua      : '))
        e = float(input('Toleransi kesalahan: '))

        if f(x0, equation) * f(x1, equation) > 0.0:
            print('Nilai tebakan awal tidak menghasilkan akar')
            print('Coba lagi dengan nilai tebakan yang berbeda')
        else:
            step = 1
            print('\n\n*** Biseksi ***')
            condition = True
            while condition:
                x2 = (x0 + x1) / 2
                fx2 = f(x2, equation)
                if fx2 is not None:
                    print('Iterasi-%d, x2 = %0.6f dan f(x2) = %0.6f' % (step, x2, fx2))
                    if f(x0, equation) * fx2 < 0:
                        x1 = x2
                    else:
                        x0 = x2
                    step += 1
                    condition = abs(fx2) > e
                else:
                    break

            if fx2 is not None:
                print('\nAkar yang dihasilkan adalah: %0.8f' % x2)

        choice = input("Lakukan perhitungan Biseksi lagi? (y/n): ")
        if choice.lower() != "y":
            break  # Keluar dari perulangan jika pengguna tidak ingin melanjutkan

# Fungsi Regula Falsi
def regulaFalsi():
    while True:  # Perulangan tambahan untuk menu Regula Falsi
        # definisi fungsi
        fx = input("Masukkan f(x): ")
        def f(x):
            return eval(fx)

        # masukkan tebakan awal dan toleransi kesalahan
        x0 = float(input('Tebakan pertama: '))
        x1 = float(input('Tebakan kedua: '))
        e = float(input('Toleransi kesalahan: '))

        # periksa kebenaran tebakan awal dan posisi falsi
        if f(x0) * f(x1) > 0.0:
            print('Nilai tebakan awal tidak menghasilkan akar')
            print('Coba lagi dengan nilai tebakan yang berbeda')
        else:
            step = 1
            print('\n\n*** Regula Falsi ***')
            condition = True
            while condition:
                x2 = x0 - (x1 - x0) * f(x0) / (f(x1) - f(x0))
                print('Iterasi-%d, x2 = %0.6f dan f(x2) = %0.6f' % (step, x2, f(x2)))

                if f(x0) * f(x2) < 0:
                    x1 = x2
                else:
                    x0 = x2

                step = step + 1
                condition = abs(f(x2)) > e

            print('\nAkar yang dihasilkan adalah: %0.8f' % x2)

        choice = input("Lakukan perhitungan Regula Falsi lagi? (y/n): ")
        if choice.lower() != "y":
            break  # Keluar dari perulangan jika pengguna tidak ingin melanjutkan

# Fungsi Newton-Raphson
def newton_raphson():
    while True:  # Perulangan untuk memungkinkan perhitungan berulang
        # Input fungsi dari pengguna
        equation = input("Masukkan fungsi f(x) contoh (4*x**3 - 15*x**2 + 17*x - 6): ")
        func = lambda x: eval(equation)

        # Define d_func as the derivative of func
        d_func = lambda x: (func(x + 1e-5) - func(x)) / 1e-5

        # Input tebakan awal (x)
        x = float(input("Masukkan tebakan awal (x0): "))

        # Input toleransi dan maksimum iterasi
        tolerance = float(input("Masukkan toleransi: "))
        max_iterations = int(input("Masukkan jumlah maksimum iterasi: "))

        if d_func(x) == 0:
            print("Newton-Raphson gagal dijalankan.")
            return None
        else:
            iterations = 1
            while abs(func(x) / d_func(x)) >= tolerance and iterations <= max_iterations:
                current_iteration_print = f"Iterasi : {iterations}, {x}"
                if func(x) == 0:
                    print(current_iteration_print + f"Solusi ditemukan : {x}")
                    return x

                x = x - func(x) / d_func(x)
                if d_func(x) == 0:
                    print("Newton-Raphson gagal dijalankan.")
                    return None
                current_iteration_print += f", {x}"
                iterations = iterations + 1

                print(current_iteration_print)

            print("\nJumlah iterasi : ", iterations)
            print("Hasil akhir : ", x)

        choice = input("Lakukan perhitungan Newton-Raphson lagi? (y/n): ")
        if choice.lower() != "y":
            break 

# Fungsi untuk dekomposisi LU
def dekomposisiLU():
    while True:  # Perulangan untuk memungkinkan perhitungan berulang
        # Input matriks A
        n = int(input("Masukkan jumlah baris dan kolom matriks A: "))

        A = []
        print("Masukkan elemen matriks A:")
        for i in range(n):
            row = []
            for j in range(n):
                row.append(float(input(f"A[{i+1}][{j+1}]: ")))
            A.append(row)

        A = np.array(A)

        # Input vektor b
        print("Masukkan vektor b:")
        b = []
        for i in range(n):
            b.append(float(input(f"b[{i+1}]: ")))

        b = np.array(b)

        # Initialize L and U
        L = np.zeros((n, n))
        U = np.zeros((n, n))

        for i in range(n):
            # Upper matrix
            for j in range(i, n):
                sum = 0
                for k in range(i):
                    sum += L[i][k] * U[k][j]
                U[i][j] = A[i][j] - sum

            # Lower matrix
            L[i][i] = 1
            for j in range(i + 1, n):
                sum = 0
                for k in range(i):
                    sum += L[j][k] * U[k][i]
                L[j][i] = (A[j][i] - sum) / U[i][i]

            # Print the current iteration
            print(f"\nIterasi {i + 1}:\n")
            print("Matriks L:")
            print(L)
            print("\nMatriks U:")
            print(U)

        # Print the final result
        print("\nMatriks L:")
        print(L)
        print("\nMatriks U:")
        print(U)

        # Solve the system using function solve_LU
        x = solve_LU(A, b)

        # Print the solution
        print("\nSolusi:")
        for i in range(len(x)):
            print(f'x_{i+1} : {x[i]}')

        choice = input("Lakukan perhitungan Dekomposisi LU lagi? (y/n): ")
        if choice.lower() != "y":
            break

def secant():
    while True:
        # Input fungsi dari pengguna
        equation = input("Masukkan fungsi f(x) (misal: x**3-x-1): ")
        f = lambda x: eval(equation)

        # Input tebakan awal (x0 dan x1)
        x0 = float(input("Masukkan tebakan awal (x0): "))
        x1 = float(input("Masukkan tebakan awal (x1): "))

        # Input toleransi, maksimum iterasi, dan batasan iterasi
        e = float(input("Masukkan toleransi (e) (misal: 0.000)1: "))
        N = int(input("Masukkan jumlah maksimum iterasi (N): "))

        print("\nIterasi  | x1          | f(x1)")
        print("-" * 30)

        for i in range(0, N):
            print("%d   | %f    | %f" % (i, x1, f(x1)))
            if abs(f(x1)) < e:
                print("\nSolusi ditemukan: x_root = %f" % x1)
                print("Dalam %d iterasi" % i)
                break

            x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
            x0 = x1
            x1 = x2

        print("\nMaksimum iterasi tercapai.")
        print("Hasil akhir: x_root = %f" % x1)

        choice = input("Lakukan perhitungan Metode Secant lagi? (y/n): ")
        if choice.lower() != "y":
            break

def gauss_jacobi(A, b, x0, tol, max_iter):
    n = len(A)
    x = x0.copy()
    x_new = np.zeros_like(x0)

    for iteration in range(max_iter):
        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * x[j]
            x_new[i] = (b[i] - sigma) / A[i][i]

        if np.allclose(x, x_new, rtol=tol):
            return x_new, iteration + 1

        x = x_new.copy()

        # Menampilkan tahap iterasi
        print(f"Iterasi {iteration + 1}:")
        print("x =", x)

    return x_new, max_iter

def menu_gauss_jacobi():
    while True:
        # Input matriks A
        n = int(input("Masukkan jumlah baris dan kolom matriks A: "))

        A = []
        print("Masukkan elemen matriks A:")
        for i in range(n):
            row = []
            for j in range(n):
                row.append(float(input(f"A[{i+1}][{j+1}]: ")))
            A.append(row)

        A = np.array(A)

        # Input vektor b
        print("Masukkan vektor b:")
        b = []
        for i in range(n):
            b.append(float(input(f"b[{i+1}]: ")))

        b = np.array(b)

        # Input tebakan awal (x0)
        x0 = np.zeros(len(b))

        # Input toleransi dan jumlah maksimum iterasi
        tolerance = float(input("Masukkan toleransi: "))
        max_iterations = int(input("Masukkan jumlah maksimum iterasi: "))

        result, iterations = gauss_jacobi(A, b, x0, tolerance, max_iterations)

        if iterations == max_iterations:
            print("Metode tidak konvergen setelah", max_iterations, "iterasi.")
        else:
            print("Solusi yang ditemukan setelah", iterations, "iterasi:")
            print(result)

        choice = input("Lakukan perhitungan Gauss-Jacobi lagi? (y/n): ")
        if choice.lower() != "y":
            break
# Menu utama
while True:
    print("╔══════════════════════════════════╗")
    print("║          METODE NUMERIK          ║")
    print("╚══════════════════════════════════╝")
    print("╔═════╦════════════════════════════╗")
    print("║ No. ║             Menu           ║")
    print("╠═════╬════════════════════════════╣")
    print("║  1  ║  Gauss Jordan              ║")
    print("║  2  ║  Biseksi                   ║")
    print("║  3  ║  Regula Falsi              ║")
    print("║  4  ║  Dekomposisi LU            ║")
    print("║  5  ║  Newton-Raphson            ║")
    print("║  6  ║  Secant                    ║")
    print("║  7  ║  Gauss-Jacobi              ║")  
    print("║  0  ║  Keluar                    ║")
    print("╚═════╩════════════════════════════╝")
    menu = input("Pilih menu: ")
    print("════════════════════════════════════")

    if menu == "1":
        gaussJordan()
    elif menu == "2":
        biseksi()
    elif menu == "3":
        regulaFalsi()
    elif menu == "4":
        dekomposisiLU()
    elif menu == "5":
        newton_raphson()
    elif menu == "6":
        secant()
    elif menu == "7":
        menu_gauss_jacobi()  # Memanggil fungsi menu Gauss-Jacobi
    elif menu == "0":
        print("Terima kasih. Keluar dari program.")
        break
    else:
        print("Pilihan menu tidak valid.")