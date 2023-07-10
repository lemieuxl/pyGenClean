"""Pipeline tree structure."""


from typing import KeysView, List, Optional

from ..report.summaries import Summary


class QCNode:
    """A node from the tree."""
    def __init__(self, name: str, *, module_name: Optional[str] = None,
                 parent: Optional[str] = None, bfile: Optional[str] = None,
                 summary: Optional[Summary] = None):
        self.name = name
        self.parent = parent
        self.data_from = set()
        self.bfile = bfile
        self.summary = summary
        self.module_name = module_name

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

    def add_summary(self, summary: Summary):
        """Add summary to the node."""
        self.summary = summary


class Tree:
    """The tree (which contain nodes)."""
    def __init__(self):
        self.tree = {}
        self.node_order = []

    def has_node(self, node: QCNode) -> bool:
        """Checks if the tree contains a node."""
        return node.name in self.tree

    def add_node(self, node: QCNode):
        """Add a node (with specified parent)."""
        # Adding the node to the tree
        assert node.name not in self.tree
        self.tree[node.name] = node
        self.node_order.append(node.name)

    def get_from_node_to_root(self, node_name: str) -> List[QCNode]:
        """Get from a node to the root of the tree."""
        nodes = []
        node = self.tree[node_name]
        while node.parent:
            nodes.append(node)
            node = self.tree[node.parent]
        nodes.append(node)

        return nodes

    def get_nodes(self) -> KeysView[str]:
        """Return all nodes from the tree."""
        return self.tree.keys()

    def get_node(self, node_name: str) -> QCNode:
        """Get the specific node."""
        return self.tree[node_name]

    def get_node_order(self) -> List[str]:
        """Get the node order"""
        return self.node_order
