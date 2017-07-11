/**
 * 
 */
package ncgop.cplex.tsl.ntu.sg;

import cplex.tsl.ntu.sg.Utility;

/**
 * @author XueYX
 *
 */
public class ConstantMatrix {

	static Double[][] V;
	static Double[][] v;
	public static void initialize(Double[][] y_up, int No)
	{
		V = new Double[No-1][No];
		for(int i =0; i< No-1;i++)
		{
			Double[] v_i = Utility.ArraySubtraction(y_up[y_up.length-1], y_up[i]);
			double norm_2 = Utility.ArrayNorm2(v_i);
			if(norm_2!=0) 
			{
				v_i = Utility.ArrayDivision(v_i,norm_2 );
			}
			V[i] = v_i;
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	public static Double[][] normalizeUtopia(Double[][] y_up, int No)
	{
		/**
		 * 
		 for i=1:No
         	if i == 1
            	v = y_up(No,:)-y_up(i,:);
        	elseif i==No
            	continue;    
        	else
            	v = [v;y_up(No,:)-y_up(i,:)];
        	end
         end
          % v: the normal vector of utopia plane
		 */
	    v = new Double[No-1][No];
		for(int i=1; i<= No; i++ )
		{
			Double[] v_i = null;
			if(i==1)
			{
				v_i =Utility.ArraySubtraction(y_up[No-1], y_up[i-1]);
			}
			else if(i==No)
			{
				continue;
			}
			else 
			{
				v_i = Utility.ArraySubtraction(y_up[No-1], y_up[i-1]);
			}
			v[i-1] = v_i;
		}
		return v;
	}
}
