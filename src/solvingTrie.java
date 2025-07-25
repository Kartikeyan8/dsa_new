import java.util.HashMap;
import java.util.Map;

public class solvingTrie {
    public class TrieNode {
        Map<Character, TrieNode> children = new HashMap<>();
        boolean isEndOfWord = false;
        TrieNode()
        {
            this.isEndOfWord = false;
        }
    }
    class Trie {
        TrieNode root;

        public Trie() {
            root = new TrieNode();
        }

        public void insert(String word) {
            TrieNode currentNode = root;
            for(Character ch : word.toCharArray()) {
                currentNode = currentNode.children.computeIfAbsent(ch,a->new TrieNode());
            }
            currentNode.isEndOfWord= true;
        }

        public boolean search(String word) {
            TrieNode currentNode = root;
            for(Character ch : word.toCharArray()) {
                currentNode = currentNode.children.get(ch);
                if(currentNode == null) {
                    return false;
                }
            }
            return currentNode.isEndOfWord;
        }

        public boolean startsWith(String prefix) {
           TrieNode node = root;
           for(Character ch : prefix.toCharArray())
           {
               node= node.children.get(ch);
               if(node == null) {
                   return false;
               }
           }
           return true;
        }
    }
}
