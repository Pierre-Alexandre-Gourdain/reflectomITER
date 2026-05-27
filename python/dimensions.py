import sympy as sp


def strong_simplify(expr):
    expr = sp.simplify(expr)
    expr = sp.together(expr)
    expr = sp.cancel(expr)
    expr = sp.factor(expr)
    expr = sp.simplify(expr)
    return expr


def parse_equation(eq_str, syms):
    lhs, rhs = eq_str.split("=", 1)
    lhs_expr = sp.sympify(lhs.strip(), locals=syms)
    rhs_expr = sp.sympify(rhs.strip(), locals=syms)
    return sp.Eq(lhs_expr, rhs_expr)


def substitute_all(equations, sol_dict):
    out = []
    for eq in equations:
        new_eq = eq.subs(sol_dict)

        if new_eq is sp.S.true or new_eq is sp.S.false:
            out.append(new_eq)
            continue

        if isinstance(new_eq, sp.Equality):
            lhs = strong_simplify(new_eq.lhs)
            rhs = strong_simplify(new_eq.rhs)
            diff = strong_simplify(lhs - rhs)

            if diff == 0:
                out.append(sp.S.true)
            else:
                out.append(sp.Eq(lhs, rhs))
        else:
            out.append(new_eq)

    return out


def solve_sequentially(equations, unknowns, free_vars, constants, solve_order):
    all_names = list(dict.fromkeys(unknowns + constants))
    syms = {name: sp.symbols(name, positive=True) for name in all_names}
    eqs = [parse_equation(eq, syms) if isinstance(eq, str) else eq for eq in equations]

    free_syms = {syms[name] for name in free_vars}
    solve_syms = [syms[name] for name in unknowns if name not in free_vars]
    solve_order_syms = [syms[name] for name in solve_order]

    for s in solve_order_syms:
        if s not in solve_syms:
            raise ValueError(f"{s} is in solve_order but is not a solved variable.")

    sol_dict = {}
    remaining_eqs = eqs[:]

    for target in solve_order_syms:
        found = False
        remaining_eqs = substitute_all(remaining_eqs, sol_dict)

        for eq in remaining_eqs:
            if eq is sp.S.true:
                continue
            if eq is sp.S.false:
                raise RuntimeError("System became inconsistent: encountered False after substitution.")
            if not isinstance(eq, sp.Equality):
                continue

            if not (eq.lhs.has(target) or eq.rhs.has(target)):
                continue

            try:
                sol_list = sp.solve(eq, target, dict=True)
            except Exception:
                sol_list = []

            if sol_list:
                expr = strong_simplify(sol_list[0][target])
                sol_dict[target] = expr
                found = True
                break

        if not found:
            raise RuntimeError(f"Could not solve explicitly for {target}")

    sol_dict = {k: strong_simplify(v.subs(sol_dict)) for k, v in sol_dict.items()}
    remaining_eqs = substitute_all(remaining_eqs, sol_dict)

    return syms, sol_dict, remaining_eqs


def print_solution_table(syms, sol_dict, unknowns, free_vars):
    free_set = {syms[v] for v in free_vars}

    print("Solution:")
    for name in unknowns:
        s = syms[name]
        if s in free_set:
            print(f"  {name} = free parameter")
        elif s in sol_dict:
            print(f"  {name} = {sp.sstr(sol_dict[s])}")
        else:
            print(f"  {name} = not solved")

def clean_physical(expr):
    expr = sp.simplify(expr)

    # Remove Abs for positive quantities
    expr = expr.replace(sp.Abs, lambda x: x)

    # Remove Heaviside (assume argument > 0 → Heaviside = 1)
    expr = expr.replace(sp.Heaviside, lambda x: 1)

    # Final simplification
    expr = sp.simplify(expr)

    return expr

from sympy.printing.cxx import cxxcode

def print_solution_defines(syms, sol_dict, unknowns, free_vars):
    free_set = {syms[v] for v in free_vars}

    print("// Free parameters (must be defined elsewhere)")
    for name in unknowns:
        if syms[name] in free_set:
            print(f"// #define {name} ...")

    print("\n// Derived constants")
    for name in unknowns:
        s = syms[name]
        if s in free_set:
            continue
        if s in sol_dict:
            expr = cxxcode(sol_dict[s])
            print(f"#define {name} ({expr})")
        else:
            print(f"// {name} not solved")

    
unknowns = ["c_0", "H_0", "B_0", "E_0", "t_0", "L_0", "T_0", "J_0", "rho_0", "p_0", "nu_0", "n_0", "m_0", "q_0"]
free_vars = ["n_0", "B_0", "L_0"]
constants = ["mu_0", "eps_0"]

eqs = [
    "E_0 / t_0 = H_0 / eps_0 / L_0",
    "E_0 / t_0 = J_0 / eps_0",
    "H_0 = B_0 / mu_0",
    "H_0 / t_0 = E_0 / mu_0 / L_0",
    "J_0 / L_0 = rho_0 / t_0",
    "L_0 / t_0 = c_0",
    "p_0 / t_0 = nu_0 * p_0",
    "p_0 / t_0 = q_0 * E_0",
    "J_0 = n_0 * q_0 * p_0 / m_0",
    "q_0 * T_0 = m_0 / eps_0 / mu_0",
    "B_0 / (m_0 * n_0 * mu_0) ** 0.5 = L_0 / t_0"
]

solve_order = ["t_0", "c_0", "E_0", "q_0", "H_0", "rho_0", "J_0", "nu_0", "p_0", "m_0", "T_0"]

syms, sol_dict, remaining_eqs = solve_sequentially(
    eqs, unknowns, free_vars, constants, solve_order
)

sol_dict = {k: clean_physical(v.subs(sol_dict)) for k, v in sol_dict.items()}

# print_solution_table(syms, sol_dict, unknowns, free_vars)
print_solution_defines(syms, sol_dict, unknowns, free_vars)

print("\nRemaining equations after substitution:")
for eq in remaining_eqs:
    print(" ", sp.sstr(eq))
    