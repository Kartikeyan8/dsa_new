import com.sun.source.tree.Tree;

import java.util.*;

public class Main {
    public static void swap(int[] nums, int i, int j) {
        int temp = nums[i];
        nums[i] = nums[j];
        nums[j] = temp;
    }

    public static <T> void view_list(int[] nums) {
        for (int i = 0; i < nums.length; i++) {
            System.out.print(nums[i] + " ");
        }
    }

    public static void permutation_helper(int[] nums, int start, ArrayList<ArrayList<Integer>> list) {
        if (start == nums.length) {
            ArrayList<Integer> temp = new ArrayList<>();
            for (Integer n : nums)
                temp.add(n);
            list.add(temp);
        }
        for (int i = start; i < nums.length; i++) {
            swap(nums, start, i);
            permutation_helper(nums, start + 1, list);
            swap(nums, start, i);
        }
    }

    public static ArrayList<ArrayList<Integer>> permutation() {
        ArrayList<ArrayList<Integer>> list = new ArrayList<>();

        int[] nums = {1, 2, 3};
        permutation_helper(nums, 0, list);
        return list;
    }

    public static void combinationSumHelper(int[] candidates, int target, List<Integer> current, List<List<Integer>> result, int start) {
        if (start == candidates.length) {
            if (target == 0) {
                result.add(new ArrayList<>(current));
                return;
            }
            return;
        }
        if (candidates[start] <= target) {
            current.add(candidates[start]);
            combinationSumHelper(candidates, target - candidates[start], current, result, start);
            current.remove(current.size() - 1);
        }
        combinationSumHelper(candidates, target, current, result, start + 1);
    }

    public static void combination_sum() {
        int[] candidates = {2, 3, 6, 7};
        int target = 7;

        Arrays.sort(candidates);

        List<List<Integer>> result = new ArrayList<>();
        combinationSumHelper(candidates, target, new ArrayList<>(), result, 0);
        System.out.println(result);
    }

    public static void sort_colors() {
        int[] nums = {2, 0, 2, 1, 1, 0};
        //0,0,2,1,1,2
        int low = 0, mid = 0, high = nums.length - 1;
        while (mid <= high) {
            if (nums[mid] == 2) {
                swap(nums, mid, high);
                high -= 1;
            } else if (nums[mid] == 0) {
                swap(nums, low, mid);
                low += 1;
                mid += 1;
            } else {
                mid += 1;
            }
        }
    }

    public static void majority_elements() {
        int[] v = {2, 1, 2, 1, 1};
        int current = v[0];
        int current_count = 1;
        for (int i = 1; i < v.length; i++) {

            if (v[i] == current) {
                current_count += 1;
            } else {
                current_count -= 1;
                if (current_count == 0) {
                    current = v[i];
                    current_count = 1;

                }
            }
        }
        System.out.println(current);
    }

    public static void prod_ex_self() {
        int[] v = {1, 2, 3, 4};
        //1,1,2,6
        int[] result = new int[v.length];
        Arrays.fill(result, 1);
        result[0] = 1;
        int left = 1;
        for (int i = 1; i < v.length; i++) {
            left *= v[i - 1];
            result[i] *= left;
        }
        int right = 1;
        for (int i = v.length - 2; i >= 0; i--) {
            right = right * v[i + 1];
            result[i] *= right;
        }

    }

    public static int find_the_dup() {
        int[] v = {1, 3, 4, 2, 2};
        Arrays.sort(v);
        for (int i = 1; i < v.length; i++) {
            if (v[i] == v[i - 1]) {
                return v[i];
            }
        }
        return -1;
    }

    public static void single_number() {
        int[] v = {4, 1, 2, 1, 2};
        int result = 0;
        for (int i = 0; i < v.length; i++) {
            result ^= v[i];
        }
        System.out.println(result);
    }

    public static void maxi_sub() {
        int[] v = {5, 4, -1, 7, 8};
        int current_max = v[0];
        int max_so_far = v[0];
        for (int i = 1; i < v.length; i++) {
            current_max = Math.max(v[i], current_max + v[i]);
            max_so_far = Math.max(max_so_far, current_max);
        }
        System.out.println(max_so_far);
    }

    public static void three_sum() {
        int[] v = {-1, 0, 1, 2, -1, -4};
        Arrays.sort(v);

        ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
        for (int i = 0; i < v.length - 2; i++) {
            if (i > 0 && v[i] == v[i - 1])
                continue;
            int left = i + 1, right = v.length - 1;
            while (left < right) {
                int current_sum = v[i] + v[left] + v[right];
                if (current_sum == 0) {
                    ans.add(new ArrayList<>(Arrays.asList(v[i], v[left], v[right])));
                    while (left < right && v[left] == v[left + 1]) {
                        left++;
                    }
                    while (right > left && v[right] == v[right - 1]) {
                        right--;
                    }
                    left += 1;
                    right -= 1;
                } else if (current_sum < 0) {
                    left += 1;
                } else {
                    right -= 1;
                }
            }


        }
        System.out.println(ans.toString());
    }

    public static int[] two_sum() {
        int[] v = {3, 2, 4};
        int target = 6;
        int[] ans = new int[2];
        Map<Integer, Integer> m = new HashMap<>();
        for (int i = 0; i < v.length; i++) {

            if (m.containsKey(target - v[i])) {
                ans[0] = m.get(target - v[i]);
                ans[1] = i;
                return ans;
            }
            m.put(v[i], i);
        }
        return ans;
    }

    public static void group_ana() {
        String[] s = {"eat", "tea", "tan", "ate", "nat", "bat"};
        Map<String, List<String>> map = new HashMap<>();
        for (String word : s) {
            char[] ch = word.toCharArray();
            Arrays.sort(ch);
            String sortedWord = new String(ch);
            if (!map.containsKey(sortedWord)) {
                map.put(sortedWord, new ArrayList<>());

            }
            map.get(sortedWord).add(word);
        }
        List<List<String>> ans = new ArrayList<>();
        System.out.println(map.values());
        for (List<String> groups : map.values()) {
            ans.add(groups);
        }

//        List<List<String>> ans = new ArrayList<>();
//        Map<String,List<String>> map = Arrays.stream(s).collect(Collectors.groupingBy(
//                p -> {
//                    char [] ch = p.toCharArray();
//                    Arrays.sort(ch);
//                    return new String(ch);
//                }
//        ));
//        for(Map.Entry<String,List<String>> entry : map.entrySet())
//        {
//            ans.add(entry.getValue());
//        }

    }

    public static int longest_cons_seq() {
        int[] v = {0, 3, 7, 2, 5, 8, 4, 6, 0, 1};
        Arrays.sort(v);
        int ans = 1;
        int res = 1;
        for (int i = 1; i < v.length; i++) {
            if (v[i] == v[i - 1] + 1) {
                ans += 1;
            } else if (v[i] != v[i - 1]) {
                ans = 1;
            }
            res = Math.max(res, ans);
        }
        return res;
    }

    public static void daily_temp() {                 //0,1, 2, 3, 4, 5, 6, 7
        int[] temp = {73, 74, 75, 71, 69, 72, 76, 73};
        Stack<Integer> st = new Stack<>();
        int[] ans = new int[temp.length];
        Arrays.fill(ans, 0);
        for (int i = temp.length - 1; i >= 0; i--) {
            while (!st.isEmpty() && temp[st.peek()] <= temp[i])
                st.pop();
            if (!st.isEmpty())
                ans[i] = st.peek() - i;
            st.add(i);
        }
        System.out.println(Arrays.toString(ans));
    }

    public static boolean valid_p() {
        String s = "({[]})";
        Stack<Character> st = new Stack<>();

        for (Character ch : s.toCharArray()) {
            if (ch == '(' || ch == '{' || ch == '[') {
                st.push(ch);
            } else {
                if (ch == ')' && !st.isEmpty()) {
                    Character top = st.pop();
                    if (top != '(')
                        return false;
                } else if (ch == '}' && !st.isEmpty()) {
                    Character top = st.pop();
                    if (top != '{')
                        return false;
                } else if (ch == ']' && !st.isEmpty()) {
                    Character top = st.pop();
                    if (top != '[')
                        return false;
                } else {
                    return false;
                }
            }
        }
        return st.isEmpty();
    }

    public static void trapping_rain_water() {
        int[] height = {4, 2, 0, 3, 2, 5};
        int[] left = new int[height.length];
        int[] right = new int[height.length];
        left[0] = height[0];
        right[height.length - 1] = height[height.length - 1];
        for (int i = 1; i < height.length; i++) {
            left[i] = Math.max(left[i - 1], height[i]);
        }
        for (int i = height.length - 2; i >= 0; i--) {
            right[i] = Math.max(right[i + 1], height[i]);
        }
        int total_water = 0;
        for (int i = 1; i < height.length - 1; i++) {
            int mini = Math.min(left[i], right[i]);
            if (mini > height[i]) {
                total_water = total_water + mini - height[i];
            }

        }
        System.out.println("Total water trapped: " + total_water);
    }

    public static void move_zeros() {
        int[] v = {0, 1, 0, 3, 0, 2, 12};
        //          12,1,0,3,0,2,0
        int left = 0, right = v.length - 1;
        for (int i = 0; i < v.length; i++) {
            if (v[i] != 0) {

                swap(v, left, i);
                left += 1;
            }
        }
        view_list(v);

    }

    public static void reverse_array(int[] v) {
        int left = 0, right = v.length - 1;
        while (left < right) {
            int temp = v[left];
            v[left] = v[right];
            v[right] = temp;
            left += 1;
            right -= 1;
        }
    }

    public static void rotate_image() {
        int[][] matrix = {
                {1, 2, 3},
                {4, 5, 6},
                {7, 8, 9}
        };
        for (int i = 0; i < matrix.length; i++) {
            for (int j = i + 1; j < matrix[0].length; j++) {
                int temp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = temp;
            }
        }
        for (int i = 0; i < matrix.length; i++) {
            reverse_array(matrix[i]);
        }
        System.out.println(Arrays.deepToString(matrix));
    }

    public static void set_zero() {
        int[][] matrix = {
                {1, 1, 1, 1},
                {1, 0, 1, 1},
                {1, 1, 1, 1}
        };
        int[] row = new int[matrix.length];
        int[] col = new int[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (matrix[i][j] == 0) {
                    row[i] = 1;
                    col[j] = 1;
                }
            }
        }
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (row[i] == 1 || col[j] == 1) {
                    matrix[i][j] = 0;
                }
            }

        }
    }

    public static void house_r() {
        int[] v = {2, 7, 9, 3, 1};
        int[] dp = new int[v.length];
        dp[0] = v[0];
        dp[1] = Math.max(v[0], v[1]);
        for (int i = 2; i < v.length; i++) {
            dp[i] = Math.max(dp[i - 1], v[i] + dp[i - 2]);
        }

    }

    public static void stairs() {
        int n = 3;
        int[] dp = new int[n];
        dp[0] = 1;
        dp[1] = 2;
        for (int i = 2; i < n; i++) {
            dp[i] = dp[i - 1] + dp[i - 2];
        }
        System.out.println(dp[n - 1]);
    }

    public static boolean jump_game() {
        int[] v = {2, 3, 1, 1, 4};
        int max_reach = 0;
        for (int i = 0; i < v.length; i++) {
            if (i > max_reach) {
                return false;
            }
            max_reach = Math.max(i + v[i], max_reach);
        }
        return true;
    }

    public static void stock1() {
        int[] prices = {7, 6, 4, 3, 1};
        int profit = 0;
        int current_min = prices[0];
        for (int i = 1; i < prices.length; i++) {
            if (prices[i] < current_min) {
                current_min = prices[i];
            } else {
                profit = Math.max(profit, prices[i] - current_min);
            }
        }
        System.out.println("Max profit: " + profit);
    }

    public static void lis() {
        int[] v = {10, 9, 2, 5, 3, 7, 101, 18};
        int[] dp = new int[v.length];
        Arrays.fill(dp, 1);
        dp[0] = 1;
        int max = 1;
        for (int i = 1; i < v.length; i++) {
            for (int j = 0; j < i; j++) {
                if (v[i] > v[j]) {
                    dp[i] = Math.max(1 + dp[j], dp[i]);
                }
                max = Math.max(max, dp[i]);
            }
        }

    }

    public static int unique_p() {
        int m = 3, n = 7;
        int[][] dp = new int[m][n];
        Arrays.fill(dp[0], 1);
        for (int i = 1; i < m; i++) {
            dp[i][0] = 1;
        }
        for (int i = 1; i < m; i++) {
            for (int j = 1; j < n; j++) {
                dp[i][j] = dp[i - 1][j] + dp[i][j - 1];
            }
        }
        return dp[m][n];
    }

    public static boolean binary_search(int[] v, int target) {
        int left = 0, right = v.length - 1;
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (v[mid] == target) {
                return true;
            } else if (v[mid] > target) {
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        return false;
    }

    public static boolean search() {
        int target = 5;
        int[][] matrix = {{1, 3, 5, 7},
                {10, 11, 16, 20},
                {23, 30, 34, 60}};
        int rows = matrix.length;
        int cols = matrix[0].length;
        for (int j = 0; j < rows; j++) {
            if (matrix[j][cols - 1] >= target) {
                System.out.println("called");
                return binary_search(matrix[j], target);
            }
        }
        return false;
    }

    public static void length_of_longest_sub() {
        String s = "abcnsdbcbb";
        int left = 0;
        int max_length = 0;
        Set<Character> set = new HashSet<>();
        for (int i = 0; i < s.length(); i++) {
            if (set.contains(s.charAt(i))) {
                while (set.contains(s.charAt(i))) {
                    set.remove(s.charAt(left));
                    left += 1;

                }
            }
            max_length = Math.max(max_length, i - left + 1);
            set.add(s.charAt(i));
        }
        System.out.println("Length of longest substring without repeating characters: " + max_length);
    }

    public static int[] sliding_win_max() {
        int[] v = {1, 3, -1, -3, 5, 3, 6, 7};
        int k = 3;
        int left = 0, right = 0;
        int[] ans = new int[v.length - k + 1];
        PriorityQueue<Integer> q = new PriorityQueue<>(Collections.reverseOrder());
        while (right < v.length) {
            q.add(v[right]);
            if (right - left + 1 == k) {
                ans[left] = q.peek();
                q.remove(v[left]);
                left += 1;
            }
            right += 1;
        }
        System.out.println(Arrays.toString(ans));
        return ans;
    }

    public static void num_of_island_helper(int[][] grid, int i, int j) {
        if (i < 0 || i >= grid.length || j < 0 || j >= grid[0].length || grid[i][j] != 1) {
            return;
        }
        grid[i][j] = 4;
        num_of_island_helper(grid, i + 1, j);
        num_of_island_helper(grid, i - 1, j);
        num_of_island_helper(grid, i, j + 1);
        num_of_island_helper(grid, i, j - 1);


    }

    public static void num_of_island() {
        int[][] grid = {
                {1, 1, 0, 0, 0},
                {1, 1, 0, 1, 1},
                {0, 0, 0, 1, 0},
                {0, 1, 1, 0, 0},
                {1, 0, 0, 0, 1}
        };
        int count = 0;

        for (int i = 0; i < grid[0].length; i++) {
            for (int j = 0; j < grid.length; j++) {
                if (grid[i][j] == 1) {
                    num_of_island_helper(grid, i, j);
                    count += 1;
                }
            }
        }

        System.out.println("Number of islands: " + count);
    }

    public static void most_water() {
        int[] height = {1, 8, 6, 2, 5, 4, 8, 3, 7};
        int left = 0, right = height.length - 1;
        int ans = 0;

        while (left < right) {
            int length = right - left;
            int breadth = Math.min(height[left], height[right]);
            ans = Math.max(ans, length * breadth);
            if (height[left] < height[right]) {
                left += 1;
            } else {
                right -= 1;
            }
        }
        System.out.println("Maximum area of water that can be contained: " + ans);
    }

    public static void sub_sum_k() {
        int[] v = {1, 2, 3};
        int k = 3;
        Map<Integer, Integer> m = new HashMap<>();
        int count = 0;
        int current_sum = 0;
        for (int i = 0; i < v.length; i++) {
            current_sum += v[i];
            if (current_sum == k) {
                count += 1;
            }
            if (m.containsKey(current_sum - k)) {
                count += m.get(current_sum - k);
            }
            m.put(current_sum, m.getOrDefault(current_sum, 0) + 1);
        }
        System.out.println("Number of subarrays with sum equal to " + k + ": " + count);
    }

    public static void subsets_helper(int[] v, int start, List<Integer> path, List<List<Integer>> ans) {
        ans.add(new ArrayList<>(path));
        for (int i = start; i < v.length; i++) {
            path.add(v[i]);
            subsets_helper(v, i + 1, path, ans);
            path.remove(path.size() - 1);
        }
    }

    public static List<List<Integer>> subsets() {
        int[] v = {1, 2, 3};
        List<List<Integer>> ans = new ArrayList<>();
        List<Integer> path = new ArrayList<>();
        subsets_helper(v, 0, path, ans);
        return ans;
    }

    public static boolean isFreshOrange(int x, int y, int[][] grid) {
        if (x < 0 || x >= grid.length || y < 0 || y >= grid[0].length || grid[x][y] != 1) {
            return false;
        }
        return true;
    }

    public static int rotting_oranges() {
        int[][] grid = {
                {2, 1, 1},
                {1, 1, 0},
                {0, 1, 1}
        };

        int n = grid.length;
        int m = grid[0].length;

        int[][] directions = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};

        int rot = 0;
        Queue<Pair<Integer, Integer>> q = new LinkedList<>();

        // Count fresh oranges and collect rotten ones
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (grid[i][j] == 1) {
                    rot++;
                }
                if (grid[i][j] == 2) {
                    q.add(new Pair<>(i, j));  // ✅ Use generics here
                }
            }
        }

        int elapsed_time = 0;
        int count = 0;

        while (!q.isEmpty()) {
            int size_of_rotten_oranges = q.size();

            for (int i = 0; i < size_of_rotten_oranges; i++) {
                Pair<Integer, Integer> current = q.poll();
                int dx = current.first;
                int dy = current.second;

                for (int[] dir : directions) {
                    int new_x = dx + dir[0];
                    int new_y = dy + dir[1];

                    if (isFreshOrange(new_x, new_y, grid)) {
                        grid[new_x][new_y] = 2;
                        q.add(new Pair<>(new_x, new_y));
                        count++;
                    }
                }
            }

            if (!q.isEmpty()) {
                elapsed_time++;
            }
        }

        return rot != count ? -1 : elapsed_time;
    }


    public static boolean word_search_helper(char[][] board, String word, int index, int i, int j) {
        if (i < 0 || i >= board.length || j < 0 || j >= board[0].length || board[i][j] != word.charAt(index)) {
            return false;
        }
        if (index == word.length() - 1) {
            return true;
        }
        char temp = board[i][j];
        board[i][j] = '*';

        boolean found = word_search_helper(board, word, index + 1, i + 1, j) ||
                word_search_helper(board, word, index + 1, i - 1, j) ||
                word_search_helper(board, word, index + 1, i, j + 1) ||
                word_search_helper(board, word, index + 1, i, j - 1);
        board[i][j] = temp;
        return found;
    }

    public static boolean word_search(char[][] board) {
        int m = board.length;
        String word = "ABCCED";
        int n = board[0].length;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (board[i][j] == word.charAt(0)) {
                    if (word_search_helper(board, word, 0, i, j)) {
                        System.out.println("Word found in the board");
                        return true;
                    }
                }
            }
        }
        return false;
    }
//        public static int min_path_sum_helper(int[][] grid, int i, int j, int m, int n, int total_sum) {
//            if()
//        }

        public static int min_path_sum() {
            int [][] grid = {
                    {1, 3, 1},
                    {1, 5, 1},
                    {4, 2, 1}
            };
            int m = grid.length;
            int n = grid[0].length;
            int [] [] dp = new int[m][n];
            Arrays.fill(dp[0], Integer.MAX_VALUE);
            for(int i = 1; i < m; i++) {
                grid[i][0] =grid[i][0] + grid[i-1][0];
            }
            for(int i = 1;i<n;i++)
            {
                grid[0][i] += grid[0][i-1];
            }
            for(int i =1;i<m;i++)
            {
                for(int j = 1;j<n;j++)
                {
                    grid[i][j] = grid[i][j] + Math.min(grid[i-1][j] , grid[i][j-1]);
                }
            }
            return grid[m-1][n-1];
        }
    public static int edit_distance_helper(String s,String t,int m,int n,int [][] dp)
    {
        if(m == 0)
        {
            return n;
        }
        if(n==0)
        {
            return m;
        }
        if(dp[m-1][n-1] != 0)
        {
            return dp[m-1][n-1];
        }
        if(s.charAt(m-1) == t.charAt(n-1))
        {
            return dp[m-1][n-1] = edit_distance_helper(s,t,m-1,n-1,dp);
        }
        return dp[m-1][n-1] = 1 + Math.min(Math.min(edit_distance_helper(s,t,m-1,n,dp),edit_distance_helper(s,t,m,n-1,dp)),edit_distance_helper(s,t,m-1,n-1,dp));
    }
    public static int edit_distance()
    {
        String s = "horse";
        String t = "ros";
        int m = s.length();
        int n = t.length();
        int [][] dp = new int[m+1][n+1];
        int a = edit_distance_helper(s,t,m,n,dp);
        System.out.println("Edit distance between \"" + s + "\" and \"" + t + "\" is: " + a);
        return a;
    }
    public static void max_prod_sub()
    {
        int [] v = {2, 3, -2, 4};
        int n = v.length;
        int lmax = 1;
        int rmax = 1;
        int ans = Integer.MIN_VALUE;
        for(int i = 0 ;i<n;i++)
        {
           if(lmax == 0)
           {
               lmax = 1;
           }
           if(rmax == 0)
           {
               rmax=1;
           }

           lmax = lmax * v[i];
           rmax = rmax * v[n-i-1];
            ans = Math.max(ans, Math.max(lmax, rmax));
        }
        System.out.println("Maximum product subarray is: " + ans);
    }
    public static class TreeNode{
        int val;
        TreeNode left;
        TreeNode right;

        TreeNode(int x) {
            val = x;
            left = null;
            right = null;
        }

    }
    public static void inorder(TreeNode root,List<Integer>ans)
    {
        if(root == null)
        {
            return ;
        }
        inorder(root.left,ans);
        ans.add(root.val);
        inorder(root.right,ans);
    }
    public static boolean symmetric_helper(TreeNode p,TreeNode q)
    {
        if(p == null && q == null)
        {
            return true;
        }
        if(p == null || q == null)
        {
            return false;
        }
        if(p.val != q.val)
        {
            return false;
        }
        return symmetric_helper(p.left,q.right) && symmetric_helper(p.right,q.left);
    }
    public static boolean symmetric(TreeNode root)
    {
        if(root == null)
        {
            return true;
        }

        return symmetric_helper(root.left,root.right);
    }
    public static List<List<Integer>> level_order(TreeNode root)
    {
        if(root == null)
        {
            return new ArrayList<>();
        }
        List<List<Integer>> ans = new ArrayList<>();
        Queue<TreeNode> q = new LinkedList<>();
        q.add(root);
        while(!q.isEmpty())
        {
            int size = q.size();
            TreeNode current = q.peek();
            List<Integer> level = new ArrayList<>();
            while(size-- > 0)
            {
                current = q.poll();
                assert current != null;
                level.add(current.val);
                if(current.left != null)
                {
                    q.add(current.left);
                }
                if(current.right != null)
                {
                    q.add(current.right);
                }

            }
            ans.add(level);

        }
        return ans;
    }
    public static int height(TreeNode root)
    {
        if(root == null)
        {
            return 0;
        }
        return 1 + Math.max(height(root.left), height(root.right));

    }
    public static List<Integer> right_side_view(TreeNode root)
    {
        if(root == null)
        {
            return new ArrayList<>();
        }
        int current_level = -1;
        Queue<Pair<TreeNode,Integer>> q = new LinkedList<>();
        List<Integer> ans = new ArrayList<>();
        q.add(new Pair<>(root,0));
        while(!q.isEmpty())
        {

            Pair<TreeNode,Integer> current = q.poll();

            if(current_level < current.second)
            {
                ans.add(current.first.val);
                current_level = current.second;
            }

            if(current.first.right != null)
            {
                q.add(new Pair<>(current.first.right,current.second + 1));
            }
            if(current.first.left != null)
            {
                q.add(new Pair<>(current.first.left,current.second + 1));
            }
        }
        System.out.println(ans);
        return ans;
    }
    public static int diameter(TreeNode root)
    {
        if(root == null)
        {
            return 0;
        }
        int left_height = height(root.left);
        int right_height = height(root.right);
        return Math.max(left_height+right_height,Math.max(diameter(root.left),diameter(root.right)));
    }
public static boolean validate_binary_search_tree_helper(TreeNode root,long left_bound,long right_bound)
{
    if(root == null)
        return true;
    if(root.val <= left_bound || root.val >= right_bound)
        return false;
    return validate_binary_search_tree_helper(root.left,left_bound,root.val) && validate_binary_search_tree_helper(root.right,root.val,right_bound);
}
public static boolean validate_binary_search_tree(TreeNode root)
{
    if(root == null)
        return true;
    if(root.left == null && root.right == null)
        return true;
    return validate_binary_search_tree_helper(root.left,Long.MIN_VALUE, Long.MAX_VALUE) ;

}
public static TreeNode invert_tree(TreeNode root)
{
    if(root == null)
    {
        return null;
    }

    TreeNode newnode = new TreeNode(root.val);
    TreeNode nnode = newnode;
    newnode.left = invert_tree(root.right);
    newnode.right = invert_tree(root.left);
    return nnode;
}
    public static void TreeNode_helper()
    {
        TreeNode root = new TreeNode(3);
        root.left = new TreeNode(9);
        root.right = new TreeNode(20);
        root.right.left = new TreeNode(15);
        root.right.right = new TreeNode(7);
        System.out.println(invert_tree(root));
    }
    public static int [][] merge_intervals()
    {
        int [][] intervals = {{1,4},{2,6},{8,10},{15,18}};
        Arrays.sort(intervals, Comparator.comparingInt(a -> a[0]));
        List<int[]> merged = new ArrayList<>();
        int [] prev = intervals[0];
        for(int i =1 ; i<intervals.length;i++)
        {
            if(intervals[i][0] < prev[1])
            {
                prev[1] = Math.max(prev[1],intervals[i][1]);
            }
            else{
                merged.add(prev);
                prev = intervals[i];
            }

        }
        merged.add(prev);
        return merged.toArray(new int[merged.size()][]);
    }
    public static void koko_eating_banana()
    {
        int [] piles = {3,6,7,11};
        int h = 8;
        int total_bananas = (Arrays.stream(piles).sum());
        int left = 1, right = total_bananas;
        while(left < right) {
            int mid = left + (right - left) / 2;
            int hours_needed = 0;
            for(int i = 0 ;i<piles.length;i++)
            {
                System.out.println();
                hours_needed = (int) (hours_needed + Math.ceil((piles[i] * 1.0) / mid)); // Using Math.ceil to round up the division
//                hours_needed += (piles[i] + mid - 1) / mid; // Equivalent to Math.ceil(piles[i] / mid)
            }
            if (hours_needed <= h) {
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        System.out.println("Minimum eating speed: " + left);
    }
    public static void first_and_last_pos()
    {
        int target = 8 ;

        int [] v = {8,8,8,8,8,8,8,8,8,8,9};
        int left = 0, right = v.length;
        int [] ans = {-1,-1};
        while(left < right)
        {
            int mid = left+(right-left)/2;
            if(v[mid] == target)
            {
                int x = mid;
                while(mid>0 && v[mid-1] == target)
                {
                    mid -= 1;
                }
                while(x<v.length && v[x] == target)
                {
                    x += 1;
                }
                ans[0] = mid;
                ans[1] = x-1 ;
                break;
            }
            else if (v[mid] <= target)
            {
                left = mid + 1;
            }
            else
            {
                right = mid - 1;
            }
        }
        System.out.println("First and last position of target " + target + " are: " + Arrays.toString(ans));
    }
    public static void jump_game2()
    {
        int [] v = {2,3,1,1,4};
        int max_reach = 0 ;
        int cur_reach = 0 ;
        int n = v.length;
        int jump = 0;
        for(int i =0 ;i<v.length;i++)
        {
            if(cur_reach < i)
            {
                jump+=1;
                cur_reach = max_reach;
            }
            max_reach = Math.max(max_reach,i+v[i]);
        }
        System.out.println("Minimum number of jumps to reach the end: " + jump);
    }
    public static int coin_change_helper(int[] coins, int target, int n) {
       if(target == 0)
        {
            return 0;
        }
         if(n == 0 || target < 0)
          {
                return Integer.MAX_VALUE;
          }
        int include = coin_change_helper(coins, target - coins[n - 1], n);
        int exclude = coin_change_helper(coins, target, n - 1);
        if(include != Integer.MAX_VALUE)
        {
            include += 1;
        }
        return Math.min(include, exclude);


    }
    public static int coin_change() {
        int[] coins = {1, 2, 5};
        int target = 11;
        int [] dp = new int [target+1];

        Arrays.fill(dp, Integer.MAX_VALUE);
        dp[0] = 0;
        for(int i = 1;i <target+1;i++)
        {
            for(int j =0 ;j<coins.length;j++)
            {
                if(i-coins[j]>=0 && dp[i - coins[j]] != Integer.MAX_VALUE)
                {
                    dp[i]= Math.min(1 + dp[i-coins[j]],dp[i]);
                }
            }

        }
        return dp[target+1];
//        System.out.println("Coin change dp: " + Arrays.toString(dp));

//        int coin = coin_change_helper(coins,target,coins.length);
//        System.out.println("Minimum number of coins needed to make " + target + " is: " + coin);
    }
    public static void perfect_square() {
        int n = 13;
        int[] dp = new int[n + 1];
        Arrays.fill(dp, Integer.MAX_VALUE);
        dp[0] = 0;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j * j <= i; j++) {
                dp[i] = Math.min(dp[i], dp[i - j * j] + 1);
            }
        }

        System.out.println("DP Table: " + Arrays.toString(dp));
        System.out.println("Minimum number of perfect squares that sum to " + n + " is: " + dp[n]);
    }
    public static int search_insert_pos()
    {
        int [] v = {1,3,5,6};
        int target = 2;
        int left = 0, right = v.length - 1;
        while(left <= right)
        {
            if(left== right)
            {
               if(target>v[left])
               {
                   return left+1;
               }
               else{
                   return left;
               }
            }
            int mid = left + (right - left) / 2;
            if(v[mid] == target)
            {
                return mid;
            }
            else if(v[mid] < target)
            {
                left = mid + 1;
            }
            else
            {
                right = mid - 1;
            }
        }
        return left;
    }
    public static int longest_common_subseq_helper(String s,String t,int m,int n,int [][]dp)
    {
        if(m<=0 || n<=0)
        {
            return 0;
        }
        if(dp[m][n] != -1)
        {
            return dp[m][n];
        }
        if(s.charAt(m-1) == t.charAt(n-1))
        {
            return dp[m][n] = 1 + longest_common_subseq_helper(s,t,m-1,n-1,dp);
        }
        return dp[m][n] = Math.max(longest_common_subseq_helper(s,t,m-1,n,dp),longest_common_subseq_helper(s,t,m,n-1,dp));

    }

    public static void longest_common_subseq()
    {
        String s1 = "abcde";
        String s2 = "ace";
        int m = s1.length();
        int n = s2.length();

        int[][] dp = new int[m + 1][n + 1];
        for (int[] row : dp) {
            Arrays.fill(row, -1);
        }
        longest_common_subseq_helper(s1,s2,m,n,dp);
        System.out.println("Length of longest common subsequence: " + dp[m][n]);
    }

        public static void spiral_matrix() {
            int[][] matrix = {
                    {1, 2, 3, 4},
                    {5, 6, 7, 8},
                    {9, 10, 11, 12},
                    {13, 14, 15, 16}
            };

            int top = 0, bottom = matrix.length - 1;
            int left = 0, right = matrix[0].length - 1;

            List<Integer> ans = new ArrayList<>();

            while (top <= bottom && left <= right) {

                // Traverse from left to right
                for (int i = left; i <= right; i++) {
                    ans.add(matrix[top][i]);
                }
                top++;

                // Traverse from top to bottom
                for (int i = top; i <= bottom; i++) {
                    ans.add(matrix[i][right]);
                }
                right--;

                // Traverse from right to left
                if (top <= bottom) {
                    for (int i = right; i >= left; i--) {
                        ans.add(matrix[bottom][i]);
                    }
                    bottom--;
                }

                // Traverse from bottom to top
                if (left <= right) {
                    for (int i = bottom; i >= top; i--) {
                        ans.add(matrix[i][left]);
                    }
                    left++;
                }
            }

            System.out.println("Spiral order of the matrix: " + ans);
        }
        public static void next_perm()
        {
            int [] v = {1,5,1};    //1,5,4,6
            int right = v.length - 1;
            int index = v.length -2;
            while(index >= 0 && v[index] >= v[index + 1])
            {
                index -= 1;
            }
            if(index <0)
            {
                Arrays.sort(v);
                System.out.println("Next permutation: " + Arrays.toString(v));
                return ;
            }
            while(right>index && v[right]<=v[index])
            {
                right-=1;
            }
            int temp = v[index];
            v[index] = v[right];
            v[right] = temp;
            System.out.println("Current permutation: " + Arrays.toString(v));
            Arrays.sort(v,index+1,v.length);
            System.out.println("Next permutation: " + Arrays.toString(v));

        }
        public static void generate_parenthesis_helper(String current,int n,int open,int close,List<String> ans)
        {
            if(open+close == 2*n)
            {
                ans.add(current);
                return;
            }
            if(open<n)
            {
                generate_parenthesis_helper(current + "(" ,n,open+1,close,ans);
            }
            if(close < open)
            {
                generate_parenthesis_helper(current + ")" ,n,open,close+1,ans);
            }

        }
        public static void generate_parenthesis()
        {
            int n =3;
            List<String> ans = new ArrayList<>();
            int open = 0, close = 0;
            generate_parenthesis_helper("",n,open,close,ans);
            System.out.println("Generated Parentheses: " + ans);
        }
        public static void letter_combination_helper(String digits, int index, String current, Map<Character, String> digitToLetters, List<String> ans) {
             if(index == digits.length())
             {

                 ans.add(current);
                 return ;
             }

                    char digit = digits.charAt(index);
                    String letters = digitToLetters.get(digit);
            System.out.println("Current digit: " + digit + ", Letters: " + letters);

                    for(int j = 0 ;j<letters.length();j++)
                    {
                        System.out.println("Here J is "+ j + " function call with current: " + current + letters.charAt(j));
                        letter_combination_helper(digits,index+1,current + letters.charAt(j),digitToLetters,ans);
                    }
    }

public static List<String> letter_combination()
{
    String digits = "23";
    Map <Character, String> digitToLetters = new HashMap<>();
    digitToLetters.put('2', "abc");
    digitToLetters.put('3', "def");
    digitToLetters.put('4', "ghi");
    digitToLetters.put('5',"jkl");
    digitToLetters.put('6', "mno");
    digitToLetters.put('7', "pqrs");
    digitToLetters.put('8', "tuv");
    digitToLetters.put('9', "wxyz");
    List<String> ans = new ArrayList<>();
    letter_combination_helper(digits,0,"",digitToLetters,ans);
    System.out.println(ans);
    return ans;


}

    public static void main(String[] args) {
        letter_combination();


    }
}

