#include "binarySearchTree.h"
#include "basicFunctions.h"
#include "utilities.h"
#include <vector>
//Global variables
int verbose = 0;
int dept = 0;

Node* insert(Node* node, Node* key, bool* alreadyThere) {
	// If the tree is empty, return a new node
    if (node == NULL){
        return key;
	}
    // If the key is already present in the tree,
    // return the node
    if (node->key == key->key){
		*alreadyThere = true;
        return node;
	}
    // Otherwise, recur down the tree/ If the key
    // to be inserted is greater than the node's key,
    // insert it in the right subtree
    if (node->key < key->key){
        node->right = insert(node->right, key,alreadyThere);
	}
    // If the key to be inserted is smaller than
    // the node's key,insert it in the left subtree
    else {
        node->left = insert(node->left, key,alreadyThere);
	}
    // Return the (unchanged) node pointer
    return node;
}

// function to search a key in a BST
const Node* search(const Node* root, uint64_t key) {

    // Base Cases: root is null or key
    // is present at root
    if (root == NULL || root->key == key)
        return root;

    // Key is greater than root's key
    if (root->key < key)
        return search(root->right, key);

    // Key is smaller than root's key
    return search(root->left, key);
}

// A utility function to do inorder tree traversal
void inorder(const Node* root) {
    if (root != NULL) {
        inorder(root->left);
		std::cout << root->key << " ";
        inorder(root->right);
    }
}

// Puts the tree in a uint64_t vector
void treeToVec(const Node* root, std::vector<sType>* vec) {
	if (root != NULL) {
		treeToVec(root->left,vec);
		vec->push_back(root->key);
		treeToVec(root->right,vec);
	}
}
void storeInorder(Node* root, std::vector<sType>* nodes) {
    if (root == nullptr)
        return;

    // Traverse the left subtree
    storeInorder(root->left, nodes);

    // Store the node data
    nodes->push_back(root->key);

    // Traverse the right subtree
    storeInorder(root->right, nodes);
}

// Function to get the count of nodes
// in complete binary tree
uint64_t totalNodes(Node* root)
{
    if (root == NULL)
        return 0;

    int l = totalNodes(root->left);
    int r = totalNodes(root->right);

    return 1 + l + r;
}


// Note that it is not a generic inorder
// successor function. It mainly works
// when right child is not empty which is
// the case wwe need in BST delete
Node* getSuccessor(Node* curr)
{
    curr = curr->right;
    while (curr != NULL && curr->left != NULL)
        curr = curr->left;
    return curr;
}

// This function deletes a given key x from
// the give BST and returns modified root of
// the BST (if it is modified)
Node* delNode(Node* root, uint64_t x)
{

    // Base case
    if (root == NULL)
        return root;

    // If key to be searched is in a subtree
    if (root->key > x)
        root->left = delNode(root->left, x);
    else if (root->key < x)
        root->right = delNode(root->right, x);

    // If root matches with the given key
    else {

        // Cases when root has 0 children
        // or only right child
        if (root->left == NULL) {
            Node* temp = root->right;
            delete root;
            return temp;
        }

        // When root has only left child
        if (root->right == NULL) {
            Node* temp = root->left;
            delete root;
            return temp;
        }

        // When both children are present
        Node* succ = getSuccessor(root);
        root->key = succ->key;
        root->right = delNode(root->right, succ->key);
    }
    return root;
}


int height(Node *root) {
    if (root == nullptr)
        return -1;

    // compute the height of left and right subtrees
    int lHeight = height(root->left);
    int rHeight = height(root->right);

    return std::max(lHeight, rHeight) + 1;
}

// Function to build a balanced BST from a sorted array
Node* buildBalancedTree(std::vector<Node>* arr, std::vector<sType>* nodes, int start, int end) {

    // Base case
    if (start > end)
        return nullptr;

    // Get the middle element and make it the root
    int mid = (start + end) / 2;
    //Node* root = new Node(nodes->at(mid));
    int loc = arr->size();
    arr->push_back(Node(nodes->at(mid)));

    // Recursively build the left and right subtrees
    arr->at(loc).left = buildBalancedTree (arr, nodes, start, mid - 1);
    arr->at(loc).right = buildBalancedTree(arr, nodes, mid + 1, end);

    return &(arr->at(loc));
}

// Function to balance a BST
void balanceBST(std::vector<Node>* arr) {
    std::vector<sType> nodes;
    // Store the nodes in sorted order
    Node root = arr->at(0);
    storeInorder(&root, &nodes);

    arr->clear();

    // Build the balanced tree from the sorted nodes
    buildBalancedTree(arr, &nodes, 0, nodes.size() - 1);
}

void printBT(const std::string& prefix, const Node* node, bool isLeft)
{
    if( node != nullptr )
    {
        std::cout << prefix;

        std::cout << (isLeft ? "├──" : "└──" );

        // print the value of the node
        std::cout << node->key << std::endl;

        // enter the next tree level - left and right branch
        printBT( prefix + (isLeft ? "│   " : "    "), node->left, node->right!=nullptr);
        printBT( prefix + (isLeft ? "│   " : "    "), node->right, false);
    }
}

void printBT(const Node* node)
{
    printBT("", node, false);
}

