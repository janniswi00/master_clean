import argparse
import pandas as pd
import networkx as nx


## Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--strain", default="PR8", type=str)
parser.add_argument("--rep", default=0, type=int)
parser.add_argument("--hpi", default=5, type=int)
parser.add_argument("--python_version", type=str)
args = parser.parse_args()

## Variables
strain = args.strain
rep = args.rep 
hpi = args.hpi
python_version = args.python_version
base_dir = "/data/influenza-genome-packaging/results/preprocessed"

## Load MST and registration data
regs = pd.read_csv(f"edges_and_registration_py{python_version}.csv").query("strain==@strain & rep==@rep & hpi==@hpi")
edges = [(source, target, weight) for source, target, weight in zip(list(regs["source"]), list(regs["target"]), list(regs["weight"]))]

# Identify connected MSTs and root for each MST
G = nx.Graph()
G.add_weighted_edges_from(edges)

for mst in nx.connected_components(G):
    sub = G.subgraph(mst)

    node_weight_sum = {
        n: sum(d.get("weight", 0) for _, _, d in sub.edges(n, data=True))
        for n in sub.nodes
    }

    root = list(node_weight_sum.keys())[list(node_weight_sum.values()).index(max(node_weight_sum.values()))]
    tree = nx.bfs_tree(sub, root)

    branches = {
        child: list(nx.descendants(tree, child)) + [child]
        for child in tree.successors(root)
    }

    # print(regs)
    # print("root: ", root)

    for branch in branches.items():
        # print("branch: ", branch)
        nodes = branch[1][::-1]
        prev_node = root
        tx = 0
        ty = 0

        for node in nodes:
            # print("node: ", node)
            # print("prev_node: ", prev_node)

            source_target_df1 = regs.query("source==@prev_node & target==@node")
            source_target_df2 = regs.query("source==@node & target==@prev_node")

            if len(source_target_df1) == 1:
                tx += source_target_df1["tx"].values[0]
                ty += source_target_df1["ty"].values[0]

            elif len(source_target_df2) == 1:
                tx += (source_target_df2["tx"].values[0] * -1)
                ty += (source_target_df2["ty"].values[0] * -1)

            # print("df1: ", source_target_df1)
            # print("df2: ", source_target_df2)
            print("Registration: ", tx, ty)

            prev_node = node