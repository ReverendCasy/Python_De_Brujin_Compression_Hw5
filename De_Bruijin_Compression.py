from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
from graphviz import Digraph
import argparse
import pydot
import random

class Vertex:

    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = defaultdict()
        self.out_edges = defaultdict()

    def increase_coverage(self):
        self.coverage += 1



class Edge:

    def __init__(self, k1, k2):
        self.k1 = k1
        self.k2 = k2
        self.seq = k1 + k2[-1]
        self.n = 2
        self.coverage = 0

    def calc_coverage(self, c1, c2):
        self.coverage = (c1 + c2) / 2

    def coverage_increment(self):
        self.coverage += 1

    def merge_next(self, next):
        self.seq +=next.seq[k:]
        self.coverage = (self.coverage * len(self.seq) + next.coverage * len(next.seq)) \
                   / (len(self.seq) + len(next.seq))

    def merge_prev(self, prev):
        self.seq = prev.seq + self.seq[k:]
        self.coverage = (self.coverage * len(self.seq) + prev.coverage * len(prev.seq)) \
                   / (len(self.seq) + len(prev.seq))

class Graph:

    def __init__(self, k):
        self.vertices = defaultdict()
        self.k = k

    def add_read(self, read):
        if len(read) < self.k:
            return

        kmer = read[:k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)

        for next_kmer_indx in range(1, len(read) - k + 1, 1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx + k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)

            new_edge = Edge(kmer, next_kmer)
            self.vertices[next_kmer].in_edges[kmer] = [new_edge]
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]
            kmer = next_kmer

    def calc_init_edge_coverage(self):

        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(
                    self.vertices[current_vertex].coverage, self.vertices[next_vertex].coverage)

    def merge_vertices(self, vertex, prev, next):
        self.vertices[prev].out_edges[next] = [Edge(prev, next)]
        self.vertices[next].in_edges[prev] = [Edge(prev, next)]
        self.vertices[prev].out_edges[next][0].merge_next(self.vertices[next].in_edges[prev][0])
        self.vertices[next].in_edges[prev][0].merge_prev(self.vertices[prev].out_edges[next][0])

        del self.vertices[vertex]
        del self.vertices[prev].out_edges[vertex]
        del self.vertices[next].in_edges[vertex]


    def compress(self):
        blacklist = []
        for vertex in self.vertices.keys():
            if (len(list(self.vertices[vertex].out_edges.keys()))) == 1 and (len(list(self.vertices[vertex].in_edges.keys())) == 1):
                blacklist.append(vertex)
        for vertex in blacklist:
            if vertex in self.vertices.keys() and len(self.vertices) > 2:
                prev, next = list(self.vertices[vertex].in_edges.keys())[0], list(self.vertices[vertex].out_edges.keys())[0]
                self.merge_vertices(vertex, prev, next)

    def cut(self, baseline):
        self.to_cut = [vertex for vertex in self.vertices.keys() if ((len(self.vertices[vertex].out_edges.keys()) == 0) and (len(self.vertices[vertex].in_edges.keys()) == 1))
                       or ((len(self.vertices[vertex].out_edges.keys()) == 1) and (len(self.vertices[vertex].in_edges.keys()) == 0))
                       or ((len(self.vertices[vertex].out_edges.keys()) == 0) and (len(self.vertices[vertex].in_edges.keys()) == 0))]
        for vertex in self.to_cut:
            if self.vertices[vertex].coverage <= baseline:
                if len(self.vertices[vertex].out_edges.keys()) == 1:
                    self.vertices[list(self.vertices[vertex].out_edges.keys())[0]].in_edges.pop(vertex)
                if len(self.vertices[vertex].in_edges.keys()) == 1:
                    self.vertices[list(self.vertices[vertex].in_edges.keys())[0]].out_edges.pop(vertex)

                del self.vertices[vertex]


    def get_contigs(self, dest):
        i = 1
        contigs = []
        for vertex in self.vertices.keys():
            for out_vertex in self.vertices[vertex].out_edges.keys():
                contigs.append(SeqRecord(Seq(self.vertices[vertex].out_edges[out_vertex][0].seq), id="contig {}".format(i)))
                i += 1
        print(contigs)
        SeqIO.write(contigs, dest, "fasta")

    def visualize(self,graph):
        vis = Digraph(comment='De Bruijin genome assemble graph')

        if graph == 'full':
            for k, v in self.vertices.items():

                vis.node(k, label=f'{k}')
                for kk, vv in v.out_edges.items():
                    vis.edge(k, kk, label=f'{vv[0].seq}')
        else:
            for k, v in self.vertices.items():

                vis.node(k, label='cov={v.coverage}')
                for kk, vv in v.out_edges.items():
                       vis.edge(k, kk, label='cov={cov} len={len}'.format(cov=vv[0].edge_coverage, len=len(vv[0].seq)))

        vis.view()
        vis.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='graph_viualization')
    parser.add_argument('-i', help='Name/path to your input fasta file', metavar='File',
                        type=str, required=True)
    parser.add_argument('-k', help='k-mer length', default=3, type=int)
    parser.add_argument('-t', help='full/not', default='full', type=str)
    parser.add_argument('-f', help='forward oriented assembly', dest='feature', action='store_true')
    parser.add_argument('-r', help='reverse complement oriented assembly', dest='feature', action='store_false')
    parser.add_argument('-c', help='Perform graph compression', action='store_true')
    parser.add_argument('-b', help='Cutoff baseline value', default=100, type=int)
    parser.add_argument('-o', help='Set path to contig output', type=str, required=False)
    parser.set_defaults(feature=True)

    args = parser.parse_args()
    i, k, t , c, b, o,  feature = args.i, args.k, args.t, args.c, args.b, args.o, args.feature
    my_graph = Graph(k)

    with open(i, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if feature == True:
                read = str(record.seq)
            else:
                read = str(record.reverse_complement().seq)
            my_graph.add_read(read)

    my_graph.calc_init_edge_coverage()
    if args.c:
        while len(my_graph.vertices) > b:
            my_graph.compress()
            my_graph.cut(2)
        if args.o:
            my_graph.get_contigs(o)

    my_graph.visualize(t)
