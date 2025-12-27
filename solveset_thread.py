import sympy
from sympy.abc import x
import multiprocessing
import time
import sys


def run_solveset(equation, symbol, domain, queue):
    """Target function to run solveset and put the result in a queue."""
    try:
        # The solveset function call is inside the separate process
        solution = sympy.solveset(equation, symbol, domain=domain)
        queue.put(solution)
    except Exception as e:
        queue.put(e)  # In case of an unexpected error


def solve_with_timeout(equation, symbol, domain=sympy.S.Reals, timeout=10):
    """Runs solveset with a timeout in a separate process."""
    # Use a Queue to get the result back from the process
    queue = multiprocessing.Queue()
    # Create a process to run the function
    process = multiprocessing.Process(
        target=run_solveset, args=(equation, symbol, domain, queue))

    # Start the process
    process.start()

    # Wait for the process to complete with a timeout
    process.join(timeout)

    # Check if the process is still alive
    if process.is_alive():
        print(f"Timeout of {timeout} seconds reached. Terminating process.")
        process.terminate()  # Forcefully terminate the process
        process.join()  # Clean up the terminated process
        return "TIMEOUT"
    else:
        # Process finished, get the result from the queue
        result = queue.get()
        if isinstance(result, Exception):
            raise result
        return result

# --- Example Usage ---


if __name__ == "__main__":
    # Example of a difficult equation that might take a long time or hang
    # (Note: SymPy behavior on specific complex equations can vary)
    equation_hard = sympy.sin(x) - sympy.sqrt(x)  # A non-trivial example

    print("Attempting to solve the equation with a 5-second timeout...")
    try:
        solution_hard = solve_with_timeout(equation_hard, x, timeout=5)
        if solution_hard == "TIMEOUT":
            print("Solving operation was interrupted due to timeout.")
        else:
            print(f"Solution: {solution_hard}")
    except Exception as e:
        print(f"An error occurred: {e}")

    print("\nAttempting to solve a simple equation...")
    # A simple equation that solves quickly
    equation_simple = x**2 - 4
    try:
        solution_simple = solve_with_timeout(
            equation_simple, x, domain=sympy.S.Reals, timeout=5)
        print(f"Solution: {solution_simple}")
    except Exception as e:
        print(f"An error occurred: {e}")
