from ete2 import Tree#,  faces

# enter the path of your newick tree here 
MyTree = "./core.tree"

def convert_to_ultrametric(root, tree_length, strategy="balanced"):
    ''' Converts a tree to ultrametric topology (all leaves have the
    same distance two root).'''

    # pre-calculate how many splits remain under each node
    node2max_depth = {}
    for node in t.traverse("postorder"):
        if not node.is_leaf():
            max_depth = max([node2max_depth[c] for c in node.children]) + 1
            node2max_depth[node] = max_depth
        else:
            node2max_depth[node] = 1
    node2dist = {root: 0.0}
    tree_length = float(tree_length)
    step = tree_length / node2max_depth[t]
    for node in t.iter_descendants("preorder"):
        if strategy == "balanced":
            node.dist = (tree_length - node2dist[node.up]) / node2max_depth[node]
            node2dist[node] =  node.dist + node2dist[node.up]
        elif strategy == "fixed":
            if not node.is_leaf():
                node.dist = step
            else:
                node.dist = tree_length - ((node2dist[node.up]) * step)
            node2dist[node] = node2dist[node.up] + 1
        node.dist = node.dist
                        
def ultrametric_layout(node):
    # node balls consume space in the tree picture, so partitions with
    # many splits are not well aligned with partitions having less
    # splits. To solve it,  I set node sphere size to 0
    node.img_style["size"] = 0
    if node.is_leaf():
        faces.add_face_to_node(nameFace, node, 0)
        
if __name__ == "__main__":
    t = Tree(MyTree)
    # Convert tree to a ultrametric topology in which distance from
    # leaf to root is always 100. Two strategies are available:
    # balanced or fixed
    convert_to_ultrametric(t, 1.0, "balanced")
    '''
    # Print distances from all leaves to root. Due to precision issues
    # with the float type.  Branch lengths may show differences at
    # high precision levels, that's way I round to 6 decimal
    # positions.
    print "distance from all leaves to root:", \
        set([round(l.get_distance(t), 6)for l in t.iter_leaves()])
    nameFace = faces.AttrFace("name")
    '''
    t.write(format=1, outfile="core_ultrametric.tree")
