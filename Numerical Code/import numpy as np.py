import numpy as np

def print_matrix(matrix):
    """Print matrix in a readable format"""
    rows, cols = matrix.shape
    for i in range(rows):
        row_str = "| "
        for j in range(cols):
            # Format to show fractions clearly
            if abs(matrix[i, j]) < 1e-10:  # Handle floating point precision
                row_str += "0.00 "
            else:
                row_str += f"{matrix[i, j]:.2f} "
        print(row_str + "|")
    print()


def gauss_elimination(A, b):
    """
    Solves Ax = b using Gauss Elimination
    
    Args:
        A: Coefficient matrix
        b: Right hand side vector
        
    Returns:
        x: Solution vector
    """
    # Create augmented matrix [A|b]
    n = len(A)
    augmented = np.column_stack((A, b))
    
    print("Initial augmented matrix:")
    print_matrix(augmented)
    
    # Forward elimination
    for i in range(n):
        # Find max pivot in current column (partial pivoting)
        max_row = i + np.argmax(abs(augmented[i:, i]))
        if max_row != i:
            augmented[[i, max_row]] = augmented[[max_row, i]]
            print(f"Swap row {i+1} with row {max_row+1}:")
            print_matrix(augmented)
        
        # Skip if pivot is zero (singular matrix)
        if abs(augmented[i, i]) < 1e-10:
            continue
            
        # Scale the pivot row to make pivot = 1
        pivot = augmented[i, i]
        augmented[i] = augmented[i] / pivot
        print(f"Scale row {i+1} by dividing by {pivot:.2f}:")
        print_matrix(augmented)
        
        # Eliminate entries below pivot
        for j in range(i + 1, n):
            factor = augmented[j, i]
            augmented[j] = augmented[j] - factor * augmented[i]
            if abs(factor) > 1e-10:  # Only print if we're actually changing something
                print(f"Subtract {factor:.2f} × row {i+1} from row {j+1}:")
                print_matrix(augmented)
    
    # Back substitution
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = augmented[i, -1]
        for j in range(i + 1, n):
            x[i] -= augmented[i, j] * x[j]
    
    print("Solution:", x)
    return x


def gauss_jordan_elimination(A, b):
    """
    Solves Ax = b using Gauss-Jordan Elimination
    
    Args:
        A: Coefficient matrix
        b: Right hand side vector
        
    Returns:
        x: Solution vector
    """
    # Create augmented matrix [A|b]
    n = len(A)
    augmented = np.column_stack((A, b))
    
    print("Initial augmented matrix:")
    print_matrix(augmented)
    
    # Forward elimination
    for i in range(n):
        # Find max pivot in current column (partial pivoting)
        max_row = i + np.argmax(abs(augmented[i:, i]))
        if max_row != i:
            augmented[[i, max_row]] = augmented[[max_row, i]]
            print(f"Swap row {i+1} with row {max_row+1}:")
            print_matrix(augmented)
        
        # Skip if pivot is zero (singular matrix)
        if abs(augmented[i, i]) < 1e-10:
            continue
            
        # Scale the pivot row to make pivot = 1
        pivot = augmented[i, i]
        augmented[i] = augmented[i] / pivot
        print(f"Scale row {i+1} by dividing by {pivot:.2f}:")
        print_matrix(augmented)
        
        # Eliminate entries below and above pivot
        for j in range(n):
            if j != i:
                factor = augmented[j, i]
                augmented[j] = augmented[j] - factor * augmented[i]
                if abs(factor) > 1e-10:  # Only print if we're actually changing something
                    print(f"Subtract {factor:.2f} × row {i+1} from row {j+1}:")
                    print_matrix(augmented)
    
    # Extract solution
    x = augmented[:, -1]
    
    print("Solution:", x)
    return x


# Example usage
if __name__ == "__main__":
    # Example system:
    # 2x + y - z = 8
    # -3x - y + 2z = -11
    # -2x + y + 2z = -3
    
    A = np.array([[2, 1, -1], 
                  [-3, -1, 2], 
                  [-2, 1, 2]], dtype=float)
    b = np.array([8, -11, -3], dtype=float)
    
    print("=== Gauss Elimination Example ===")
    print("Solving the system:")
    print("2x + y - z = 8")
    print("-3x - y + 2z = -11")
    print("-2x + y + 2z = -3")
    print()
    
    gauss_elimination(A.copy(), b.copy())
    
    print("\n=== Gauss-Jordan Elimination Example ===")
    print("Solving the same system:")
    print("2x + y - z = 8")
    print("-3x - y + 2z = -11")
    print("-2x + y + 2z = -3")
    print()
    
    gauss_jordan_elimination(A.copy(), b.copy())