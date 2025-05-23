import tkinter as tk
from tkinter import messagebox, scrolledtext

def center_window(win, width=500, height=400):
    screen_width = win.winfo_screenwidth()
    screen_height = win.winfo_screenheight()
    x = (screen_width // 2) - (width // 2)
    y = (screen_height // 2) - (height // 2)
    win.geometry(f"{width}x{height}+{x}+{y}")
#---------------------------------------------------- Gauss Elimination ---------------------------------------------------
#---------------------------------------------------- Gauss Elimination ---------------------------------------------------
#---------------------------------------------------- Gauss Elimination ---------------------------------------------------

def gaussian_elimination(A, b):
    n = len(A)
    steps = "Starting Gaussian Elimination:\n"
    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(A[r][i]))
        if i != max_row:
            A[i], A[max_row] = A[max_row], A[i]
            b[i], b[max_row] = b[max_row], b[i]
            steps += f"Swapped row {i+1} with row {max_row+1}\n"
        pivot = A[i][i]
        if abs(pivot) < 1e-12:
            return None, "No unique solution (pivot is zero)", steps
        steps += f"Pivot at row {i+1}: {pivot:.4f}\n"
        for j in range(i, n):
            A[i][j] /= pivot
        b[i] /= pivot
        steps += f"Normalized row {i+1}: {['{:.4f}'.format(x) for x in A[i]]}, b = {b[i]:.4f}\n"
        for k in range(i+1, n):
            factor = A[k][i]
            for j in range(i, n):
                A[k][j] -= factor * A[i][j]
            b[k] -= factor * b[i]
            steps += f"Eliminated row {k+1} using row {i+1}, factor = {factor:.4f}\n"
    x = [0] * n
    steps += "\nStarting Back Substitution:\n"
    for i in range(n-1, -1, -1):
        x[i] = b[i] - sum(A[i][j] * x[j] for j in range(i+1, n))
        steps += f"x{i+1} = {x[i]:.4f}\n"
    return x, None, steps

def gaus_elimination_solver():
    size = 3
    solver_window = tk.Toplevel(root)
    solver_window.title("Gaussian Elimination Solver")
    center_window(solver_window, 600, 450)

    entries_A = []
    tk.Label(solver_window, text="Enter Matrix A:").grid(row=0, column=0, columnspan=size)
    for i in range(size):
        row_entries = []
        for j in range(size):
            entry = tk.Entry(solver_window, width=6)
            entry.grid(row=i+1, column=j, padx=2, pady=2)
            row_entries.append(entry)
        entries_A.append(row_entries)

    tk.Label(solver_window, text="Enter Vector b:").grid(row=0, column=size+1)
    entry_b = []
    for i in range(size):
        entry = tk.Entry(solver_window, width=6)
        entry.grid(row=i+1, column=size+1, padx=2, pady=2)
        entry_b.append(entry)

    output_text = scrolledtext.ScrolledText(solver_window, width=70, height=15, font=("Courier", 10))
    output_text.grid(row=size+2, column=0, columnspan=size+2, pady=10, padx=10)

    def solve():
        try:
            A = []
            for i in range(size):
                row = []
                for j in range(size):
                    val = entries_A[i][j].get()
                    row.append(float(val))
                A.append(row)
            b = [float(entry_b[i].get()) for i in range(size)]

            import copy
            A_copy = copy.deepcopy(A)
            b_copy = copy.deepcopy(b)
            solution, err, steps = gaussian_elimination(A_copy, b_copy)
            if err:
                messagebox.showerror("Error", err, parent=solver_window)
            else:
                result_text = steps + "\nSolution:\n" + "\n".join([f"x{i+1} = {val:.4f}" for i, val in enumerate(solution)])
                output_text.delete('1.0', tk.END)
                output_text.insert(tk.END, result_text)
        except Exception:
            messagebox.showerror("Input Error", "Please enter valid numbers.", parent=solver_window)

    solve_button = tk.Button(solver_window, text="Solve", command=solve)
    solve_button.grid(row=size+1, column=0, columnspan=size+2, pady=5)

#---------------------------------------------------- Jordan  Elimination ---------------------------------------------------
def jordan_elimination(A, b):
    n = len(A)
    steps = "Starting Jordan Elimination:\n"
    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(A[r][i]))
        if i != max_row:
            A[i], A[max_row] = A[max_row], A[i]
            b[i], b[max_row] = b[max_row], b[i]
            steps += f"Swapped row {i+1} with row {max_row+1}\n"
        pivot = A[i][i]
        if abs(pivot) < 1e-12:
            return None, "No unique solution (pivot is zero)", steps
        steps += f"Pivot at row {i+1}: {pivot:.4f}\n"
        for j in range(i, n):
            A[i][j] /= pivot
        b[i] /= pivot
        steps += f"Normalized row {i+1}: {['{:.4f}'.format(x) for x in A[i]]}, b = {b[i]:.4f}\n"
        for k in range(n):
            if k != i:
                factor = A[k][i]
                for j in range(i, n):
                    A[k][j] -= factor * A[i][j]
                b[k] -= factor * b[i]
                steps += f"Eliminated row {k+1} using row {i+1}, factor = {factor:.4f}\n"
    x = [0] * n
    steps += "\nSolution:\n"
    for i in range(n):
        x[i] = b[i]
        steps += f"x{i+1} = {x[i]:.4f}\n"
    return x, None, steps
def jordan_elimination_solver():
    size = 3
    solver_window = tk.Toplevel(root)
    solver_window.title("Jordan Elimination Solver")
    center_window(solver_window, 600, 450)

    entries_A = []
    tk.Label(solver_window, text="Enter Matrix A:").grid(row=0, column=0, columnspan=size)
    for i in range(size):
        row_entries = []
        for j in range(size):
            entry = tk.Entry(solver_window, width=6)
            entry.grid(row=i+1, column=j, padx=2, pady=2)
            row_entries.append(entry)
        entries_A.append(row_entries)

    tk.Label(solver_window, text="Enter Vector b:").grid(row=0, column=size+1)
    entry_b = []
    for i in range(size):
        entry = tk.Entry(solver_window, width=6)
        entry.grid(row=i+1, column=size+1, padx=2, pady=2)
        entry_b.append(entry)

    output_text = scrolledtext.ScrolledText(solver_window, width=70, height=15, font=("Courier", 10))
    output_text.grid(row=size+2, column=0, columnspan=size+2, pady=10, padx=10)

    def solve():
        try:
            A = []
            for i in range(size):
                row = []
                for j in range(size):
                    val = entries_A[i][j].get()
                    row.append(float(val))
                A.append(row)
            b = [float(entry_b[i].get()) for i in range(size)]

            import copy
            A_copy = copy.deepcopy(A)
            b_copy = copy.deepcopy(b)
            solution, err, steps = jordan_elimination(A_copy, b_copy)
            if err:
                messagebox.showerror("Error", err, parent=solver_window)
            else:
                result_text = steps + "\nSolution:\n" + "\n".join([f"x{i+1} = {val:.4f}" for i, val in enumerate(solution)])
                output_text.delete('1.0', tk.END)
                output_text.insert(tk.END, result_text)
        except Exception:
            messagebox.showerror("Input Error", "Please enter valid numbers.", parent=solver_window)

    solve_button = tk.Button(solver_window, text="Solve", command=solve)
    solve_button.grid(row=size+1, column=0, columnspan=size+2, pady=5)



#---------------------------------------------------- Lu decomposition ---------------------------------------------------
#---------------------------------------------------- Lu decomposition ---------------------------------------------------
#---------------------------------------------------- Lu decomposition ---------------------------------------------------

def lu_decomposition(A):
    n = len(A)
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]
    steps = "Starting LU Decomposition:\n"
    for i in range(n):
        L[i][i] = 1
        for j in range(i, n):
            U[i][j] = A[i][j]
            for k in range(i):
                U[i][j] -= L[i][k] * U[k][j]
        for j in range(i+1, n):
            L[j][i] = A[j][i]
            for k in range(i):
                L[j][i] -= L[j][k] * U[k][i]
            L[j][i] /= U[i][i]
        steps += f"Step {i+1}: L = {L}, U = {U}\n"
    return L, U, steps
def lu_decomposition_solver():
    size = 3
    solver_window = tk.Toplevel(root)
    solver_window.title("LU Decomposition Solver")
    center_window(solver_window, 600, 450)

    entries_A = []
    tk.Label(solver_window, text="Enter Matrix A:").grid(row=0, column=0, columnspan=size)
    for i in range(size):
        row_entries = []
        for j in range(size):
            entry = tk.Entry(solver_window, width=6)
            entry.grid(row=i+1, column=j, padx=2, pady=2)
            row_entries.append(entry)
        entries_A.append(row_entries)

    output_text = scrolledtext.ScrolledText(solver_window, width=70, height=15, font=("Courier", 10))
    output_text.grid(row=size+1, column=0, columnspan=size+2, pady=10, padx=10)

    def solve():
        try:
            A = []
            for i in range(size):
                row = []
                for j in range(size):
                    val = entries_A[i][j].get()
                    row.append(float(val))
                A.append(row)

            import copy
            A_copy = copy.deepcopy(A)
            L, U, steps = lu_decomposition(A_copy)
            result_text = steps + "\nL:\n" + "\n".join(["\t".join([f"{val:.4f}" for val in row]) for row in L]) + \
                "\n\nU:\n" + "\n".join(["\t".join([f"{val:.4f}" for val in row]) for row in U])
            output_text.delete('1.0', tk.END)
            output_text.insert(tk.END, result_text)
        except Exception:
            messagebox.showerror("Input Error", "Please enter valid numbers.", parent=solver_window)

    solve_button = tk.Button(solver_window, text="Solve", command=solve)
    solve_button.grid(row=size+1, column=0, columnspan=size+2, pady=5)




#---------------------------------------------------- Cramer's Rule ---------------------------------------------------
#---------------------------------------------------- Cramer's Rule ---------------------------------------------------
#---------------------------------------------------- Cramer's Rule ---------------------------------------------------


def cramer_rule_full_solution(A, b):
    import copy
    n = len(A)
    det_A = determinant(A)
    steps = f"Determinant of A = {det_A:.4f}\n"
    if abs(det_A) < 1e-12:
        return None, "No unique solution (determinant is zero)", steps
    x = []
    for i in range(n):
        A_copy = copy.deepcopy(A)
        for row in range(n):
            A_copy[row][i] = b[row]
        det_i = determinant(A_copy)
        x_i = det_i / det_A
        x.append(x_i)
        steps += f"Determinant of A with column {i+1} replaced: {det_i:.4f} => x{i+1} = {x_i:.4f}\n"
    return x, None, steps

def determinant(matrix):

    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]
    det = 0
    for c in range(n):
        minor = [[matrix[r][cc] for cc in range(n) if cc != c] for r in range(1, n)]
        det += ((-1)**c) * matrix[0][c] * determinant(minor)
    return det

def cramer_solver():
    size = 3
    solver_window = tk.Toplevel(root)
    solver_window.title("Cramer's Rule Solver")
    center_window(solver_window, 600, 450)

    entries_A = []
    tk.Label(solver_window, text="Enter Matrix A:").grid(row=0, column=0, columnspan=size)
    for i in range(size):
        row_entries = []
        for j in range(size):
            entry = tk.Entry(solver_window, width=6)
            entry.grid(row=i+1, column=j, padx=2, pady=2)
            row_entries.append(entry)
        entries_A.append(row_entries)

    tk.Label(solver_window, text="Enter Vector b:").grid(row=0, column=size+1)
    entry_b = []
    for i in range(size):
        entry = tk.Entry(solver_window, width=6)
        entry.grid(row=i+1, column=size+1, padx=2, pady=2)
        entry_b.append(entry)

    output_text = scrolledtext.ScrolledText(solver_window, width=70, height=15, font=("Courier", 10))
    output_text.grid(row=size+2, column=0, columnspan=size+2, pady=10, padx=10)

    def solve():
        try:
            A = []
            for i in range(size):
                row = []
                for j in range(size):
                    val = entries_A[i][j].get()
                    row.append(float(val))
                A.append(row)
            b = [float(entry_b[i].get()) for i in range(size)]

            solution, err, steps = cramer_rule_full_solution(A, b)
            if err:
                messagebox.showerror("Error", err, parent=solver_window)
            else:
                result_text = steps + "\nSolution:\n" + "\n".join([f"x{i+1} = {val:.4f}" for i, val in enumerate(solution)])
                output_text.delete('1.0', tk.END)
                output_text.insert(tk.END, result_text)
        except Exception:
            messagebox.showerror("Input Error", "Please enter valid numbers.", parent=solver_window)

    solve_button = tk.Button(solver_window, text="Solve", command=solve)
    solve_button.grid(row=size+1, column=0, columnspan=size+2, pady=5)

#---------------------------------------------------- Main Page ---------------------------------------------------
#---------------------------------------------------- Main Page ---------------------------------------------------

root = tk.Tk()
root.title("Main Window")
center_window(root, 400, 350)

label = tk.Label(root, text="Choose a solver:")
label.pack(pady=10)

btn_gauss = tk.Button(root, text="Gaussian Elimination Solver", command=gaus_elimination_solver)
btn_gauss.pack(pady=5)

btn_cramer = tk.Button(root, text="Cramer's Rule Solver", command=cramer_solver)
btn_cramer.pack(pady=5)

btn_jordan = tk.Button(root, text="Jordan Elimination Solver", command=jordan_elimination_solver)
btn_jordan.pack(pady=5)

btn_lu = tk.Button(root, text="LU Decomposition Solver", command=lu_decomposition_solver)
btn_lu.pack(pady=5)

root.mainloop()
