import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collector;

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
    public static void main(String[] args) {
    prod_ex_self();
    }
}