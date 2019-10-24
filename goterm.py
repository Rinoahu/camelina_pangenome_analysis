#!usr/bin/env python
import networkx as nx
from networkx.algorithms.traversal.depth_first_search import dfs_tree
import sys
import os
import re
import gzip


p = re.compile(r'GO:\d{7,7}')


stanza = set(['[Term]', '[Typedef]', '[Instance]'])

# go annotation dict build
def goannot(f):
	annot_dict = {}
	goid, name = None, None
	for i in f:
		if i.startswith('id:'):
			if goid and name:
				annot_dict[goid] = name
			goid = i[:-1].split(' ')[1]

		elif i.startswith('name:'):
			name = i[6:-1]
		else:
			continue

	if goid and name:
		annot_dict[goid] = name

	return annot_dict


# the go term parse function
def oboparse(f):
	node = {}
	for i in f:
		#if i.startswith(term):
		flag = i.strip()
		#if flag in stanza:
		if flag.startswith('['):
			if node.get('stanza') == '[Term]':
				yield node
			node = {'stanza': flag}
		else:
			j = i.split(':', 1)
			if len(j) != 2:
				continue
			try:
				node[j[0]].append(j[1])
			except:
				node[j[0]] = [j[1]]

	if node.get('stanza') == '[Term]':
		yield node



# fetch the obo file for GO and save to go database
www = 'http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo'
#if not os.path.isfile('/home/zhans/tools/hx_genome_script/GOterm_KEGG/gene_ontology.1_2.obo.gz'):
if 0:
	os.system('wget -c ' + www)
	os.system('gzip --best gene_ontology.1_2.obo')
	os.system('mv gene_ontology.1_2.obo /home/zhans/tools/hx_genome_script/GOterm_KEGG/')

# find all path from source to target:
# nx.all_simple_paths(G, source=0, target=3)

# find the all node with distance less than cutoff(depath) from the source
# nx.single_source_shortest_path_length(G, source, cutoff)

graph = nx.DiGraph()
#g = nx.Graph()
#f = gzip.open('/home/zhans/tools/hx_genome_script/GOterm_KEGG/gene_ontology.1_2.obo.gz', 'r')
f = open('go.obo', 'r')
goterms = oboparse(f)
tables = {}
for goterm in goterms:
#	print 'goterm', goterm
	if 'is_obsolete' in goterm:
		continue
	goid = p.findall(goterm['id'][0])[0]
	if goid not in graph:
		graph.add_node(goid)
	tables[goid] = goterm
	# link the is_a to current term
	is_as = [[p.findall(elem)[0], goid] for elem in goterm.get('is_a', [])]
	if is_as:
		graph.add_edges_from(is_as)
	else:
		graph.add_edge('root', goid)

# print 'node', goterm
#print 'the graph', graph

# give the go id and level, return the go id at given level
def golv(goid, level = 1, goterms = None):
	if level < 1:
		level = 1
	if not goterms:
		return 'unknown'
	try:
		path = nx.shortest_path(graph, 'root', goid)
#	try:
		return path[level]
	except:
		return 'unknown'

if __name__ == '__main__':

	if len(sys.argv[1:]) < 1:
		print 'python this.py foo.txt'
		print 'foo.txt contain 2 cols seperated by tab in the following format'
		print 'locus	goid'
		raise SystemExit()

	f = open(sys.argv[1], 'r')
	outputs = {}
	visit = set()
	for i in f:
		goids = p.findall(i)
		qid = i[: -1].split('\t')[0]
		for goid in goids:
			if (qid, goid) not in visit:
				visit.add((qid, goid))
			else:
				continue

			sid = golv(goid, 2, goterms)
			if sid == 'unknown':
				continue
			namespace, name = tables[sid]['namespace'][0].strip(), tables[sid]['name'][0].strip()
			
  			output = [qid, sid, namespace, name]
#			print '\t'.join(output)
			try:
				outputs[(namespace, name)] += 1
			except:
				outputs[(namespace, name)] = 1


	keys = outputs.keys()
	keys.sort()
	for key in keys:
		print '\t'.join(key) + '\t' + str(outputs[key])

