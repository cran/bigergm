from infomap import Infomap
import numpy

def infomap_python(matrix, n_clusters):
  im = Infomap("--flow-model undirected --silent --preferred-number-of-modules "+str(n_clusters))
  
  im.add_links(matrix)
  im.run()
  res = []
  for node in im.tree:
    if node.is_leaf:
      res.append(node.module_id)
  return res


def infomap_python_directed(matrix, n_clusters):
  im = Infomap("--flow-model directed --silent --preferred-number-of-modules "+str(n_clusters))
  im.add_links(matrix)
  im.run()
  res = []
  for node in im.tree:
    if node.is_leaf:
      res.append(node.module_id)
  return res


