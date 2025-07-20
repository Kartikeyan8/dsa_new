import java.util.*;

public class Main {
    public static void swap(int[] nums, int i, int j) {
        int temp = nums[i];
        nums[i] = nums[j];
        nums[j] = temp;
    }
    public static <T> void view_list(int[] nums) {
        for(int i = 0; i < nums.length; i++) {
            System.out.print(nums[i] + " ");
        }
    }
    public static void permutation_helper(int[] nums, int start, ArrayList<ArrayList<Integer>>list) {
        if(start == nums.length)
        {
            ArrayList<Integer> temp = new ArrayList<>();
            for(Integer n : nums)
                temp.add(n);
            list.add(temp);
        }
        for(int i = start;i<nums.length;i++)
        {
            swap(nums,start,i);
            permutation_helper(nums,start+1,list);
            swap(nums,start,i);
        }
    }
    public static ArrayList<ArrayList<Integer>> permutation()
    {
        ArrayList<ArrayList<Integer>> list = new ArrayList<>();

        int [] nums = {1, 2, 3};
        permutation_helper(nums,0,list);
        return list;
    }
    public static void combinationSumHelper(int[] candidates, int target, List<Integer> current, List<List<Integer>> result, int start) {
        if(start == candidates.length)
        {
            if(target == 0)
            {
                result.add(new ArrayList<>(current));
                return;
        }
            return;
        }
        if(candidates[start]<=target)
        {
            current.add(candidates[start]);
            combinationSumHelper(candidates, target - candidates[start], current, result, start);
            current.remove(current.size() - 1);
        }
        combinationSumHelper(candidates, target, current, result, start+1);
    }
public static void combination_sum()
{
    int[] candidates = {2, 3, 6, 7};
    int target = 7;

    Arrays.sort(candidates);

    List<List<Integer>> result = new ArrayList<>();
    combinationSumHelper(candidates, target, new ArrayList<>(), result, 0);
    System.out.println(result);
}
public static void sort_colors()
{
    int[] nums = {2,0,2,1,1,0};
    //0,0,2,1,1,2
    int low = 0, mid = 0, high = nums.length - 1;
    while (mid<=high)
    {
        if(nums[mid] == 2)
        {
            swap(nums, mid, high);
            high-=1;
        }
        else if (nums[mid] == 0)
        {swap(nums, low, mid);
            low+=1;
            mid+=1;
        }
        else
        {
            mid+=1;
        }
    }
}

    public static void majority_elements() {
        int[] v = {2, 1,2,1,1};
        int current = v[0];
        int current_count = 1;
        for (int i = 1; i < v.length; i++) {

            if (v[i] == current) {
                current_count += 1;
            } else {
                current_count -= 1;
                if(current_count == 0)
                {
                    current = v[i];
                    current_count = 1;

                }
            }
        }
        System.out.println(current);
    }
    public static void prod_ex_self()
    {
        int [] v= {1,2,3,4};
        //1,1,2,6
        int [] result = new int[v.length];
        Arrays.fill(result,1);
        result[0] = 1;
        int left = 1;
        for(int i =1;i<v.length;i++)
        {
            left *= v[i-1];
            result[i] *= left;
        }
        int right = 1;
        for(int i =v.length-2;i>=0;i--)
        {
            right = right * v[i+1];
            result[i] *= right;
        }

    }
    public static int find_the_dup() {
        int[] v = {1,3,4,2,2};
        Arrays.sort(v);
        for(int i= 1;i<v.length;i++)
        {
            if(v[i] == v[i-1])
            {
                return v[i];
            }
        }
        return -1;
    }
    public static void single_number() {
        int[] v = {4,1,2,1,2};
        int result = 0;
        for(int i = 0;i<v.length;i++)
        {
            result ^= v[i];
        }
        System.out.println(result);
    }
    public static void maxi_sub()
    {
        int [] v = {5,4,-1,7,8};
        int current_max = v[0];
        int max_so_far = v[0];
        for(int i=1;i<v.length;i++)
        {
            current_max = Math.max(v[i],current_max+v[i]);
            max_so_far = Math.max(max_so_far,current_max);
        }
        System.out.println(max_so_far);
    }
    public static void three_sum()
    {
        int [] v = {-1,0,1,2,-1,-4};
        Arrays.sort(v);

        ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
        for(int i= 0;i<v.length-2;i++)
        {
            if(i>0 && v[i] == v[i-1])
                continue;
            int left = i+1,right= v.length-1;
            while(left<right)
            {
                int current_sum = v[i] + v[left] + v[right];
                if(current_sum ==0 )
                {
                    ans.add(new ArrayList<>(Arrays.asList(v[i],v[left],v[right])));
                    while(left<right && v[left] == v[left+1])
                    {
                        left++;
                    }
                    while (right>left && v[right] == v[right-1])
                    {
                        right--;
                    }
                    left+=1;
                    right-=1;
                }
                 else if (current_sum < 0)
                    {
                        left+=1;
                    }
                    else
                    {
                        right-=1;
                    }
            }


        }
        System.out.println(ans.toString());
    }
    public static int[] two_sum()
    {
        int [] v ={3,2,4};
        int target = 6;
        int[] ans = new int[2];
        Map<Integer,Integer>m =  new HashMap<>();
        for(int i =0 ;i <v.length;i++)
        {

            if(m.containsKey( target - v[i]))
            {
                ans[0] = m.get( target - v[i]);
                ans[1] = i;
                return ans;
            }
            m.put(v[i],i);
        }
        return ans;
    }
    public static void group_ana()
    {
        String [] s = {"eat", "tea", "tan", "ate", "nat", "bat"};
        Map<String,List<String>> map = new HashMap<>();
        for(String word : s)
        {
            char [] ch = word.toCharArray();
            Arrays.sort(ch);
            String sortedWord = new String (ch);
            if(!map.containsKey(sortedWord))
            {
                map.put(sortedWord,new ArrayList<>());

            }
            map.get(sortedWord).add(word);
        }
        List<List<String>> ans = new ArrayList<>();
        System.out.println(map.values());
        for(List<String> groups : map.values())
        {
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
    public static int longest_cons_seq()
    {
        int [] v = {0,3,7,2,5,8,4,6,0,1};
        Arrays.sort(v);
        int ans = 1;
        int res = 1;
        for(int i = 1;i<v.length;i++)
        {
            if(v[i] == v[i-1] + 1)
            {
                ans += 1;
            }
            else if (v[i] != v[i-1])
            {
                ans = 1;
            }
            res = Math.max(res,ans);
        }
        return res;
    }
    public static void daily_temp()
    {                 //0,1, 2, 3, 4, 5, 6, 7
        int [] temp = {73,74,75,71,69,72,76,73};
        Stack<Integer>st = new Stack<>();
        int [] ans = new int[temp.length];
        Arrays.fill(ans,0);
        for(int i = temp.length-1;i>=0;i--)
        {
            while(!st.isEmpty() && temp[st.peek()] <= temp[i])
                st.pop();
            if(!st.isEmpty())
                ans[i] = st.peek() - i;
            st.add(i);
        }
        System.out.println(Arrays.toString(ans));
    }
    public static boolean valid_p()
    {
        String s = "({[]})";
        Stack<Character> st = new Stack<>();

        for(Character ch: s.toCharArray() )
        {
            if (ch == '(' || ch == '{' || ch == '[') {
                st.push(ch);
            }
            else {
                if(ch == ')' && !st.isEmpty())
                {
                    Character top = st.pop();
                    if(top != '(')
                        return false;
                }
                else if (ch == '}' && !st.isEmpty())
                {
                    Character top = st.pop();
                    if(top != '{')
                        return false;
                }
                else if (ch == ']' && !st.isEmpty())
                {
                    Character top = st.pop();
                    if(top != '[')
                        return false;
                }
                else
                {
                    return false;
                }
            }
        }
        return st.isEmpty();
    }
    public static void trapping_rain_water()
    {
        int [] height = {4,2,0,3,2,5};
        int [] left = new int[height.length];
        int [] right = new int[height.length];
        left[0] = height[0];
        right[height.length-1] = height[height.length-1];
        for(int i=1;i<height.length;i++)
        {
            left[i] = Math.max(left[i-1],height[i]);
        }
        for(int i = height.length -2 ;i>=0;i--)
        {
            right[i] = Math.max(right[i+1],height[i]);
        }
        int total_water = 0 ;
        for(int i =1;i<height.length-1;i++)
        {
            int mini = Math.min(left[i],right[i]);
            if(mini > height[i])
            {
                total_water = total_water + mini - height[i];
            }

        }
        System.out.println("Total water trapped: " + total_water);
    }
    public static void move_zeros()
    {
        int [] v = {0,1,0,3,0,2,12};
        //          12,1,0,3,0,2,0
        int left = 0,right = v.length-1;
        for(int i =0 ;i<v.length;i++)
        {
            if(v[i] != 0)
            {

                swap(v,left,i);
                left+=1;
            }
        }
        view_list(v);

        }
        public static void reverse_array(int[] v)
        {
            int left= 0,right =v.length-1;
            while(left<right)
            {
                int temp = v[left];
                v[left] = v[right];
                v[right] = temp;
                left+=1;
                right-=1;
            }
        }
        public static void rotate_image()
        {
            int [][] matrix = {
                    {1,2,3},
                    {4,5,6},
                    {7,8,9}
            };
            for(int i =0;i<matrix.length;i++)
            {
                for(int j = i+1;j<matrix[0].length;j++)
                {
                    int temp = matrix[i][j];
                    matrix[i][j] = matrix[j][i];
                    matrix[j][i] =temp;
                }
            }
            for(int i =0 ;i<matrix.length;i++)
            {
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
                    if(row[i] == 1 || col[j] == 1)
                    {
                        matrix[i][j] = 0;
                    }
                }

            }
        }
        public static void house_r()
        {
            int [] v = {2,7,9,3,1};
            int [] dp = new int[v.length];
            dp[0] = v[0];
            dp[1] = Math.max(v[0],v[1]);
            for(int i =2;i<v.length;i++)
            {
                dp[i] = Math.max(dp[i-1],v[i] + dp[i-2]);
            }

        }
        public static void stairs()
        {
            int n = 3;
            int [] dp = new int[n];
            dp[0] =1;
            dp[1] = 2;
            for(int i =2 ;i<n;i++)
            {
                dp[i] = dp[i-1] + dp[i-2];
            }
            System.out.println(dp[n-1]);
        }
        public static boolean jump_game()
        {
            int [] v ={2,3,1,1,4};
            int max_reach = 0;
            for(int i =0 ;i<v.length;i++)
            {
                if(i>max_reach)
                {
                    return false;
                }
                max_reach = Math.max(i+v[i],max_reach);
            }
            return true;
        }
        public static void stock1()
        {
            int [] prices = {7,6,4,3,1};
            int profit = 0;
            int current_min = prices[0];
            for(int i  =1; i<prices.length;i++)
            {
                if(prices[i] < current_min)
                {
                    current_min = prices[i];
                }
                else{
                    profit = Math.max(profit, prices[i] - current_min);
                }
            }
            System.out.println("Max profit: " + profit);
        }
        public static void lis()
        {
            int [] v = {10,9,2,5,3,7,101,18};
            int [] dp = new int [v.length];
            Arrays.fill(dp,1);
            dp[0] = 1;
            int max = 1;
            for(int i =1;i<v.length;i++)
            {
                for(int j =0;j<i;j++)
                {
                    if(v[i]>v[j])
                    {
                        dp[i] = Math.max(1+dp[j],dp[i]);
                    }
                    max = Math.max(max,dp[i]);
                }
            }

        }
        public static int unique_p()
        {
            int m =3,n=7;
            int[][] dp = new int[m][n];
            Arrays.fill(dp[0],1);
            for(int i = 1; i < m; i++) {
                dp[i][0] = 1;
            }
            for(int i = 1;i<m;i++)
            {
                for(int j =1;j<n;j++)
                {
                    dp[i][j] = dp[i-1][j] + dp[i][j-1];
                }
            }
            return dp[m][n];
        }
        public static boolean binary_search(int [] v,int target)
        {
            int left = 0, right = v.length-1;
            while(left<=right)
            {
                int mid = left + (right - left) / 2;
                if(v[mid] == target)
                {
                    return true;
                }
                else if(v[mid] > target)
                {
                    right = mid - 1;
                }
                else
                {
                    left = mid + 1;
                }
            }
            return false;
        }
        public static boolean search()
        {
            int target = 5;
            int [][] matrix = {{1,3,5,7},
                               {10,11,16,20},
                               {23,30,34,60}};
            int rows = matrix.length;
            int cols = matrix[0].length;
            for(int j = 0; j<rows; j++)
            {
                if(matrix[j][cols-1] >= target)
                {
                    System.out.println("called");
                    return binary_search(matrix[j],target);
                }
            }
                return false;
        }
        public static void length_of_longest_sub()
        {
            String s = "abcnsdbcbb";
            int left = 0;
            int max_length = 0;
            Set<Character> set = new HashSet<>();
            for(int i = 0; i<s.length();i++)
            {
                if(set.contains(s.charAt(i)))
                {
                    while(set.contains(s.charAt(i)))
                    {
                        set.remove(s.charAt(left));
                        left+=1;

                    }
                }
                max_length = Math.max(max_length, i - left + 1);
                set.add(s.charAt(i));
            }
            System.out.println("Length of longest substring without repeating characters: " + max_length);
        }
        public static int [] sliding_win_max()
        {
            int [] v = {1,3,-1,-3,5,3,6,7};
            int k =3;
            int left = 0,right = 0;
            int [] ans = new int[v.length - k + 1];
            PriorityQueue<Integer> q = new PriorityQueue<>(Collections.reverseOrder());
            while(right<v.length)
            {
                q.add(v[right]);
                if(right - left + 1 == k)
                {
                    ans[left] = q.peek();
                    q.remove(v[left]);
                    left+=1;
                }
                right+=1;
            }
            System.out.println(Arrays.toString(ans));
            return ans;
        }
        public static void num_of_island_helper(int [][] grid, int i, int j)
        {
            if(i<0 || i>= grid.length || j< 0 || j>= grid[0].length || grid[i][j] != 1)
            {
                return ;
            }
            grid[i][j] = 4;
            num_of_island_helper(grid,i+1,j);
            num_of_island_helper(grid,i-1,j);
            num_of_island_helper(grid,i,j+1);
            num_of_island_helper(grid,i,j-1);


        }
        public static void num_of_island()
        {
            int [][] grid = {
                    {1,1,0,0,0},
                    {1,1,0,1,1},
                    {0,0,0,1,0},
                    {0,1,1,0,0},
                    {1,0,0,0,1}
            };
            int count = 0;

            for(int i = 0;i<grid[0].length;i++)
            {
                for(int j = 0;j<grid.length;j++)
                {
                    if(grid[i][j] == 1)
                    {
                        num_of_island_helper(grid,i,j);
                        count+=1;
                    }
                }
            }

            System.out.println("Number of islands: " + count);
        }
        public static void most_water()
        {
            int [] height = {1,8,6,2,5,4,8,3,7};
            int left =0 , right = height.length-1;
            int ans =0 ;

            while(left<right)
            {
                int length = right-left;
                int breadth = Math.min(height[left],height[right]);
                ans = Math.max(ans, length * breadth);
                if(height[left] < height[right])
                {
                    left+=1;
                }
                else
                {
                    right-=1;
                }
            }
            System.out.println("Maximum area of water that can be contained: " + ans);
        }
        public static void sub_sum_k()
        {
            int [] v = {1,2,3};
            int k = 3;
            Map<Integer,Integer> m = new HashMap<>();
            int count = 0;
            int current_sum = 0;
            for(int i = 0;i<v.length;i++)
            {
                current_sum += v[i];
                if(current_sum == k)
                {
                    count+=1;
                }
                if(m.containsKey(current_sum - k))
                {
                    count += m.get(current_sum - k);
                }
                m.put(current_sum, m.getOrDefault(current_sum, 0) + 1);
            }
            System.out.println("Number of subarrays with sum equal to " + k + ": " + count);
        }
        public static void subsets_helper(int [] v, int start, List<Integer> path, List<List<Integer>> ans)
        {
            ans.add(new ArrayList<>(path));
            for(int i = start;i<v.length;i++)
            {
                path.add(v[i]);
                subsets_helper(v,i+1,path,ans);
                path.remove(path.size()-1);
            }
        }
        public static List<List<Integer>> subsets()
        {
            int [] v = {1,2,3};
            List<List<Integer>> ans = new ArrayList<>();
            List<Integer> path = new ArrayList<>();
            subsets_helper(v,0,path,ans);
            return ans;
        }
    public static void main(String[] args) {
        subsets();
    }
}