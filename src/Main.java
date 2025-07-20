import java.util.*;
import java.util.stream.Collector;
import java.util.stream.Collectors;

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
    public static void rain()
    {
//        int [] height = {4,2,0,3,2,5};
//        int [] left = new int[height.length];
//        int [] right = new int[height.length];
//        left[0] = height[0];
//        right[height.length-1] = height[height.length-1];
//        for(int i =1;i<height.length;i++)
//        {
//            left[i] = Math.max(left[i-1],height[i]);
//        }
//        for(int i = height.length-2;i>=0;i--)
//        {
//            right[i] = Math.max(right[i+1],height[i]);
//        }
//        System.out.println("Left: " + Arrays.toString(left));
//        System.out.println("Right: " + Arrays.toString(right));
//        int ans = 0;
//        for(int i = 0;i<height.length;i++)
//        {
//            ans += Math.min(left[i],right[i]) - height[i];
//        }
//        System.out.println("Total water trapped: " + ans);
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
    public static void main(String[] args) {
        stairs();
    }
}