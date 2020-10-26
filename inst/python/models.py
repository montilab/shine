import networkx as nx
from networkx.generators.community import LFR_benchmark_graph

def LFR(seed, n, tau1, tau2, mu, average_degree, min_community, max_community, **kwargs):
    for s in range(int(seed), int(seed)+250):
        print(s)
        try:
            print("Tring seed...{0}:  ".format(s), end="")
            G = LFR_benchmark_graph(n=int(n),
                                    tau1=tau1,
                                    tau2=tau2,
                                    mu=mu,
                                    max_iters=250,
                                    average_degree=average_degree,
                                    min_community=int(min_community),
                                    max_community=int(max_community),
                                    seed=s,
                                    **kwargs)
            print("Success!")
            return(nx.adjacency_matrix(G))
        except Exception as e:
            print(e)

    raise ValueError("Maximum search iterations reached before finding a graph.")