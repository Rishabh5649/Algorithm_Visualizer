import tkinter as tk
from tkinter import messagebox, simpledialog
import matplotlib.pyplot as plt
import networkx as nx

# ---------------------- m-COLORABILITY ----------------------

def is_safe_color(node, color, colors, graph):
    for neighbor in graph[node]:
        if colors[neighbor] == color:
            return False
    return True

def color_graph(graph, m, colors, node=0):
    if node == len(graph):
        return True
    for c in range(1, m + 1):
        if is_safe_color(node, c, colors, graph):
            colors[node] = c
            if color_graph(graph, m, colors, node + 1):
                return True
            colors[node] = 0
    return False

def m_colorability_interface(root):
    win = tk.Toplevel(root)
    win.title("m-Colorability Visualizer")
    win.geometry("850x600")
    win.configure(bg="#2e2e2e")

    tk.Label(win, text="Left Click: Add Node | Right Click: Add Edge", fg="orange", bg="#2e2e2e", font=("Arial", 12)).pack(pady=5)

    canvas = tk.Canvas(win, width=600, height=500, bg="#1e1e1e")
    canvas.pack(side=tk.LEFT, padx=10, pady=10)

    control_frame = tk.Frame(win, bg="#2e2e2e")
    control_frame.pack(side=tk.RIGHT, padx=10)

    graph = {}
    nodes = []
    node_radius = 18
    selected_node = [None]

    def add_node(event):
        idx = len(nodes)
        x, y = event.x, event.y
        nodes.append((x, y))
        graph[idx] = []
        canvas.create_oval(x - node_radius, y - node_radius, x + node_radius, y + node_radius,
                           fill="#87CEFA", outline="white", width=2)
        canvas.create_text(x, y, text=str(idx), font=("Helvetica", 10, "bold"), fill="white")

    def get_node(x, y):
        for i, (nx, ny) in enumerate(nodes):
            if (x - nx) ** 2 + (y - ny) ** 2 <= node_radius ** 2:
                return i
        return None

    def add_edge(event):
        clicked = get_node(event.x, event.y)
        if clicked is None:
            return
        if selected_node[0] is None:
            selected_node[0] = clicked
            canvas.itemconfig("current", fill="yellow")
        else:
            if clicked != selected_node[0]:
                x1, y1 = nodes[selected_node[0]]
                x2, y2 = nodes[clicked]
                canvas.create_line(x1, y1, x2, y2, fill="cyan", width=2)
                if clicked not in graph[selected_node[0]]:
                    graph[selected_node[0]].append(clicked)
                if selected_node[0] not in graph[clicked]:
                    graph[clicked].append(selected_node[0])
            selected_node[0] = None

    def find_chromatic_number(graph):
        n = len(graph)
        for m in range(1, n + 1):
            colors = [0] * n
            if color_graph(graph, m, colors):
                return m, colors
        return n, list(range(1, n + 1))  # fallback (shouldn't happen with proper graphs)

    def check_m_coloring():
        if not graph:
            messagebox.showerror("Error", "Draw a graph first!")
            return

        chromatic_num, chromatic_colors = find_chromatic_number(graph)
        messagebox.showinfo("Chromatic Number", f"Minimum number of colors needed: {chromatic_num}")

        # Ask user if they want to test a custom m
        user_m = simpledialog.askinteger("Number of Colors", f"Enter value of m to test (≥ {chromatic_num}):",
                                         minvalue=1)
        if user_m is None:
            return

        colors = [0] * len(graph)
        if color_graph(graph, user_m, colors):
            messagebox.showinfo("Result", f"The graph IS {user_m}-colorable.")
            G = nx.Graph()
            for u in graph:
                for v in graph[u]:
                    G.add_edge(u, v)

            pos = nx.spring_layout(G, seed=42)
            node_colors = [colors[node] for node in G.nodes()]
            cmap = plt.cm.get_cmap('tab10', max(node_colors))
            nx.draw(G, pos, with_labels=True, node_color=node_colors, cmap=cmap, node_size=500, font_color="white",
                    edge_color="gray")
            plt.title(f"{user_m}-Coloring")
            plt.show()
        else:
            messagebox.showinfo("Result", f"The graph is NOT {user_m}-colorable.")

    def clear_all():
        canvas.delete("all")
        graph.clear()
        nodes.clear()
        selected_node[0] = None

    canvas.bind("<Button-1>", add_node)
    canvas.bind("<Button-3>", add_edge)

    tk.Button(control_frame, text="Check Coloring", command=check_m_coloring,
              font=("Arial", 12), bg="#4caf50", fg="white", width=20).pack(pady=10)
    tk.Button(control_frame, text="Clear Graph", command=clear_all,
              font=("Arial", 12), bg="#f44336", fg="white", width=20).pack(pady=10)

# ---------------------- N-QUEENS ----------------------

def n_queens_solve(n):
    board = [-1] * n
    results = []

    def is_valid(row, col):
        for i in range(row):
            if board[i] == col or abs(board[i] - col) == abs(i - row):
                return False
        return True

    def solve(row):
        if row == n:
            results.append(board[:])
            return
        for col in range(n):
            if is_valid(row, col):
                board[row] = col
                solve(row + 1)
                board[row] = -1

    solve(0)
    return results

def n_queens_interface(root):
    win = tk.Toplevel(root)
    win.title("N-Queens Visualizer")
    win.geometry("800x800")
    win.configure(bg="#2e2e2e")

    tk.Label(win, text="Enter number of queens (N ≥ 4):", font=("Arial", 14), bg="#2e2e2e", fg="white").pack(pady=10)

    input_frame = tk.Frame(win, bg="#2e2e2e")
    input_frame.pack(pady=5)
    entry = tk.Entry(input_frame, font=("Arial", 16), width=10, justify='center')
    entry.pack()
    entry.focus()

    canvas_container = tk.Frame(win, bg="#2e2e2e")
    canvas_container.pack()

    solution_buttons = []

    def on_submit():
        for widget in canvas_container.winfo_children():
            widget.destroy()

        try:
            n = int(entry.get())
            if n < 4:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid number ≥ 4")
            return

        solutions = n_queens_solve(n)
        if not solutions:
            messagebox.showinfo("Result", f"No solutions found for N = {n}")
            return

        canvas = tk.Canvas(canvas_container, width=60 * n, height=60 * n)
        canvas.pack()

        idx = [0]

        def draw_board(solution):
            canvas.delete("all")
            for i in range(n):
                for j in range(n):
                    color = "#444" if (i + j) % 2 == 0 else "#888"
                    canvas.create_rectangle(j * 60, i * 60, (j + 1) * 60, (i + 1) * 60, fill=color)
            for i in range(n):
                canvas.create_oval(solution[i]*60 + 15, i*60 + 15, solution[i]*60 + 45, i*60 + 45, fill="white")

        def next_solution():
            idx[0] = (idx[0] + 1) % len(solutions)
            draw_board(solutions[idx[0]])

        draw_board(solutions[0])

        btn = tk.Button(canvas_container, text="Next Solution", font=("Arial", 12), command=next_solution,
                        bg="#4caf50", fg="white")
        btn.pack(pady=10)
        solution_buttons.clear()
        solution_buttons.append(btn)

    tk.Button(win, text="Submit", command=on_submit, font=("Arial", 12), bg="#2196f3", fg="white").pack(pady=5)
# ---------------------- BRANCH AND BOUND - KNAPSACK ----------------------
import heapq

class Node:
    def __init__(self, level, profit, weight, bound, items_taken):
        self.level = level
        self.profit = profit
        self.weight = weight
        self.bound = bound
        self.items_taken = items_taken

    def __lt__(self, other):
        return self.bound > other.bound  # max-heap

def bound(u, n, W, items):
    if u.weight >= W:
        return 0
    profit_bound = u.profit
    j = u.level + 1
    totweight = u.weight

    while j < n and totweight + items[j][1] <= W:
        totweight += items[j][1]
        profit_bound += items[j][0]
        j += 1

    if j < n:
        profit_bound += (W - totweight) * items[j][0] / items[j][1]

    return profit_bound

def branch_and_bound_knapsack(W, weights, values):
    n = len(values)
    items = sorted(zip(values, weights), key=lambda x: x[0] / x[1], reverse=True)
    q = []
    v = Node(-1, 0, 0, 0.0, [])
    v.bound = bound(v, n, W, items)
    heapq.heappush(q, v)
    max_profit = 0
    best_items = []

    steps = []

    while q:
        v = heapq.heappop(q)
        steps.append(f"Exploring node at level {v.level}, profit = {v.profit}, weight = {v.weight}, bound = {v.bound:.2f}")
        if v.bound > max_profit:
            u = Node(0, 0, 0, 0.0, [])
            u.level = v.level + 1
            if u.level >= n:
                continue

            # Include the item
            u.weight = v.weight + items[u.level][1]
            u.profit = v.profit + items[u.level][0]
            u.items_taken = v.items_taken + [u.level]
            if u.weight <= W and u.profit > max_profit:
                max_profit = u.profit
                best_items = u.items_taken.copy()
            u.bound = bound(u, n, W, items)
            heapq.heappush(q, u)

            # Exclude the item
            u2 = Node(u.level, v.profit, v.weight, 0.0, v.items_taken.copy())
            u2.bound = bound(u2, n, W, items)
            heapq.heappush(q, u2)

    return max_profit, [i + 1 for i in best_items], steps  # Item indices start from 1
def branch_and_bound_interface(root):
    win = tk.Toplevel(root)
    win.title("Branch and Bound - Knapsack")
    win.geometry("850x600")
    win.configure(bg="#2e2e2e")

    tk.Label(win, text="Enter weights and values as comma-separated numbers", font=("Arial", 12), fg="white", bg="#2e2e2e").pack(pady=5)
    tk.Label(win, text="Example: Weights: 2,3,4 | Values: 3,4,5 | Capacity: 5", font=("Arial", 10), fg="orange", bg="#2e2e2e").pack(pady=5)

    input_frame = tk.Frame(win, bg="#2e2e2e")
    input_frame.pack(pady=10)

    tk.Label(input_frame, text="Weights:", fg="white", bg="#2e2e2e").grid(row=0, column=0, padx=5)
    weights_entry = tk.Entry(input_frame, width=30)
    weights_entry.grid(row=0, column=1, padx=5)

    tk.Label(input_frame, text="Values:", fg="white", bg="#2e2e2e").grid(row=1, column=0, padx=5)
    values_entry = tk.Entry(input_frame, width=30)
    values_entry.grid(row=1, column=1, padx=5)

    tk.Label(input_frame, text="Capacity:", fg="white", bg="#2e2e2e").grid(row=2, column=0, padx=5)
    capacity_entry = tk.Entry(input_frame, width=30)
    capacity_entry.grid(row=2, column=1, padx=5)

    output_box = tk.Text(win, height=18, width=90, bg="#1e1e1e", fg="white", font=("Courier", 10))
    output_box.pack(pady=10)

    def run_bnb():
        output_box.delete("1.0", tk.END)
        try:
            weights = list(map(int, weights_entry.get().split(",")))
            values = list(map(int, values_entry.get().split(",")))
            W = int(capacity_entry.get())
            if len(weights) != len(values) or W <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid integer lists and positive capacity.")
            return

        max_profit, items_taken, steps = branch_and_bound_knapsack(W, weights, values)
        output_box.insert(tk.END, "\n".join(steps))
        output_box.insert(tk.END, f"\n\nMaximum Profit: {max_profit}")
        output_box.insert(tk.END, f"\nItems Taken (1-based index): {items_taken}")

    tk.Button(win, text="Solve", command=run_bnb, bg="#4caf50", fg="white", font=("Arial", 12), width=20).pack()

# ---------------------- HAMILTONIAN CYCLE ----------------------

# [Unchanged from current version for brevity; same dark theme and instruction label additions can be done if desired]
def hamiltonian_cycle_interface(root):
    win = tk.Toplevel(root)
    win.title("Hamiltonian Cycle Visualizer")
    win.geometry("850x600")
    win.configure(bg="#2e2e2e")

    tk.Label(win, text="Left Click: Add Node | Right Click: Add Edge", fg="orange", bg="#2e2e2e", font=("Arial", 12)).pack(pady=5)

    canvas = tk.Canvas(win, width=600, height=500, bg="#1e1e1e")
    canvas.pack(side=tk.LEFT, padx=10, pady=10)

    control_frame = tk.Frame(win, bg="#2e2e2e")
    control_frame.pack(side=tk.RIGHT, padx=10)

    graph = {}
    nodes = []
    node_radius = 18
    selected_node = [None]

    def add_node(event):
        idx = len(nodes)
        x, y = event.x, event.y
        nodes.append((x, y))
        graph[idx] = []
        canvas.create_oval(x - node_radius, y - node_radius, x + node_radius, y + node_radius,
                           fill="#f06292", outline="white", width=2)
        canvas.create_text(x, y, text=str(idx), font=("Helvetica", 10, "bold"), fill="white")

    def get_node(x, y):
        for i, (nx, ny) in enumerate(nodes):
            if (x - nx) ** 2 + (y - ny) ** 2 <= node_radius ** 2:
                return i
        return None

    def add_edge(event):
        clicked = get_node(event.x, event.y)
        if clicked is None:
            return
        if selected_node[0] is None:
            selected_node[0] = clicked
            canvas.itemconfig("current", fill="yellow")
        else:
            if clicked != selected_node[0]:
                x1, y1 = nodes[selected_node[0]]
                x2, y2 = nodes[clicked]
                canvas.create_line(x1, y1, x2, y2, fill="orange", width=2)
                if clicked not in graph[selected_node[0]]:
                    graph[selected_node[0]].append(clicked)
                if selected_node[0] not in graph[clicked]:
                    graph[clicked].append(selected_node[0])
            selected_node[0] = None

    def is_hamiltonian_cycle(path, graph):
        if len(path) != len(graph):
            return False
        last = path[-1]
        return path[0] in graph[last]

    def hamiltonian_util(path, visited):
        if len(path) == len(graph):
            return is_hamiltonian_cycle(path, graph)
        for v in graph[path[-1]]:
            if not visited[v]:
                visited[v] = True
                path.append(v)
                if hamiltonian_util(path, visited):
                    return True
                path.pop()
                visited[v] = False
        return False

    def check_hamiltonian():
        if not graph:
            messagebox.showerror("Error", "Draw a graph first!")
            return
        n = len(graph)
        visited = [False] * n
        path = [0]
        visited[0] = True
        if hamiltonian_util(path, visited):
            messagebox.showinfo("Result", f"Hamiltonian cycle exists: {' -> '.join(map(str, path + [path[0]]))}")
        else:
            messagebox.showinfo("Result", "No Hamiltonian cycle found.")

    def clear_all():
        canvas.delete("all")
        graph.clear()
        nodes.clear()
        selected_node[0] = None

    canvas.bind("<Button-1>", add_node)
    canvas.bind("<Button-3>", add_edge)

    tk.Button(control_frame, text="Check Hamiltonian Cycle", command=check_hamiltonian,
              font=("Arial", 12), bg="#4caf50", fg="white", width=22).pack(pady=10)
    tk.Button(control_frame, text="Clear Graph", command=clear_all,
              font=("Arial", 12), bg="#f44336", fg="white", width=22).pack(pady=10)

# ---------------------- MAIN MENU ----------------------

def main_menu():
    root = tk.Tk()
    root.title("Algorithm Visualizer")
    root.geometry("700x500")
    root.configure(bg="#1e1e1e")

    tk.Label(root, text="Algorithm Visualizer", font=("Helvetica", 24, "bold"), fg="white", bg="#1e1e1e").pack(pady=20)

    problem_frame = tk.Frame(root, bg="#1e1e1e")
    problem_frame.pack(pady=10)

    def create_problem_section(name, about_text, start_func):
        frame = tk.Frame(problem_frame, bg="#333", bd=2, relief="groove", padx=10, pady=10)
        frame.pack(fill="x", pady=10, padx=40)

        tk.Label(frame, text=name, font=("Helvetica", 18, "bold"), bg="#333", fg="#00e676").pack(anchor="w")

        btn_frame = tk.Frame(frame, bg="#333")
        btn_frame.pack(anchor="w", pady=5)

        tk.Button(btn_frame, text="Start", font=("Helvetica", 12), width=12, bg="#4caf50", fg="white",
                  command=lambda: start_func(root)).pack(side="left", padx=5)

        def show_about():
            messagebox.showinfo(f"About {name}", about_text)

        tk.Button(btn_frame, text="About", font=("Helvetica", 12), width=12, bg="#2196f3", fg="white",
                  command=show_about).pack(side="left", padx=5)

    create_problem_section(
        "m-Colorability",
        "Check if a graph can be colored with a specific number of colors without any adjacent nodes sharing the same color.",
        m_colorability_interface
    )
    create_problem_section(
        "N-Queens",
        "Find the possible solutions to the N-Queens problem where no queens threaten each other.",
        n_queens_interface
    )
    create_problem_section(
        "Hamiltonian Cycle",
        "Find if a Hamiltonian cycle exists in a graph (a cycle that visits every node exactly once).",
        hamiltonian_cycle_interface
    )
    create_problem_section(
        "Branch and Bound (Knapsack)",
        "Solve the 0/1 Knapsack problem using the Branch and Bound method with upper-bound pruning.",
        branch_and_bound_interface
    )

    root.mainloop()

if __name__ == "__main__":
    main_menu()
