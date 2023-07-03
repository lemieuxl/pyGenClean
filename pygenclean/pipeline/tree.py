"""Pipeline tree structure."""


class QCNode:
    """A node from the tree."""
    def __init__(self, name: str, parent: str = None, bfile: str = None):
        self.name = name
        self.parent = parent
        self.data_from = set()
        self.bfile = bfile

    def __repr__(self) -> str:
        """String representation"""
        return (
            f"Node({self.name}, parent={self.parent}, "
            f"data_from={self.data_from})"
        )

    def change_parent(self, parent: str):
        """Change the parent."""
        self.parent = parent

    def add_data_from_node(self, node: str):
        """This current node get data from another one (except bfile)."""
        self.data_from.add(node)

    def set_bfile(self, bfile: str):
        """Set the bfile value (if any)"""
        self.bfile = bfile


class Tree:
    """The tree (which contain nodes)."""
    def __init__(self):
        self.tree = {}

    def has_node(self, node: QCNode):
        """Checks if the tree contains a node."""
        return node.name in self.tree

    def add_node(self, node: QCNode):
        """Add a node (with specified parent)."""
        # Adding the node to the tree
        self.tree[node.name] = node

    def get_from_node_to_root(self, node_name: str):
        """Get from a node to the root of the tree."""
        nodes = []
        node = self.tree[node_name]
        while node.parent:
            nodes.append(node)
            node = self.tree[node.parent]
        nodes.append(node)

        return nodes

    def get_nodes(self):
        """Return all nodes from the tree."""
        return self.tree.keys()
