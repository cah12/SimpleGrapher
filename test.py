from sympy import symbols, preorder_traversal

x, y = symbols('x y')
expr = (2**x + x*y)**(1/2)

print("Preorder Traversal:")
for subtree in preorder_traversal(expr):
    print(subtree)
