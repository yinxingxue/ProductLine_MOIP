package ncgop.cplex.tsl.ntu.sg;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import matrix.cplex.tsl.ntu.sg.src.matrix.Matrix;
import matrix.cplex.tsl.ntu.sg.src.matrix.MatrixMathematics;

import org.apache.commons.lang3.ArrayUtils;

import cplex.tsl.ntu.sg.CplexResult;
import cplex.tsl.ntu.sg.Utility;

public class SolRep_CutTry {
	/** Input*/
	Double[][] y_up;
	Double[][] f;
    int varNo;
	int n;
	Vector<LinkedHashMap<Short, Double>> ori_A;
	Vector<Double> ori_B;
	
	Vector<LinkedHashMap<Short, Double>> ori_Aeq;
	Vector<Double> ori_Beq;
	
	Vector<LinkedHashMap<Short, Double>> extra_A;
	Vector<Double> extra_B;
	
	Vector<LinkedHashMap<Short, Double>> extra_Aeq;
	Vector<Double> extra_Beq;
	
	public Double[][] getY_up() {
		return y_up;
	}

	public int getN() {
		return n;
	}

	/** Output*/
	Vector<Double[]> P;
	IloNumVar[] xVar;
	Double[] w;
	
	protected IloCplex cplex;
	protected double[] lb;
	protected double[] ub;
	
	public SolRep_CutTry(Double[][] y_up, int varNo, int n) {
		this.y_up = y_up;
		this.varNo = varNo;
		this.n = n;
		this.P = new Vector<Double[]>();
	}

	@SuppressWarnings({ "unchecked" })
	public Vector<Double[]> calculateWandV() throws Exception {
	 	int No = y_up.length;
		int Nv = this.varNo;
			
		/**  Matlab code :
		  ini = rand(1,No);
	      ini = ini/sum(ini);
  	       X0 = ini*y_up;                    % get an initial point on utopia plane
        p_k = X0;                           % put X0 into p
		 */
		
		// do the initialization of X0
		Double[] ini = Utility.generateRandomDouble(No);
		//ini[0] = 0.4970845481049563; ini[1] = 0.061224489795918366;ini[2]= 0.4380466472303207;ini[3]= 0.0036443148688046646;
		ini = Utility.ArrayDivision(ini, Utility.arraySum(ini));
		Double[] X0 = Utility.ArrayProduceMatrix(ini, Utility.MatrixTranspose(y_up));
		Double[] p = X0;
		P.add(p);
	 
	 
		/**
		 *  % Find the vectors
    	V = zeros(1,No);
    	for i=1:(No-1)
        	Vi = y_up(No,:)-y_up(i,:);
        	if norm(Vi) ~=0
            	Vi = Vi/norm(Vi);
        	end;
        	V = [V;Vi];
    	end
    	V = V(2:No,:);
    	% V: the vectors of utopia plane
		 */
		// generate constant matrix V
		ConstantMatrix.initialize(y_up, No);
		Double[][] V = ConstantMatrix.V;
		
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
		ConstantMatrix.normalizeUtopia(y_up, No);
		Double[][] v = ConstantMatrix.v;
		
		/**
		 *     
		 AAA = v(:,2:size(v,2));
    	 AAA(1,1)=5;
    	 BBB = -v(:,1);
    	 w = AAA\BBB;
    	 w = [1;w];
    	 w = w/sqrt(sum(w.^2));		
		 */
		Double[][] AAA = Utility.MatrixSubMat(v,0,1);		
		AAA[0][0] = 5.0;
		Double[] BBB = Utility.negArray(Utility.MatrixTranspose(v)[0]);
		double[][] BBB_data = new double[BBB.length][1];
		for(int i =0; i< BBB.length;i++)
		{
			BBB_data[i][0] =  BBB[i];
		}
		Matrix AAAMtr = new Matrix(Utility.toPrimateMat(AAA));
		Matrix AAAInv = MatrixMathematics.inverse(AAAMtr);
		Matrix wMtr = MatrixMathematics.multiply(AAAInv, new Matrix(BBB_data));
		double[] w_values = MatrixMathematics.transpose(wMtr).getValues()[0];
		double[] w_ext= new double[w_values.length+1]; 
		w_ext[0]= 1;
		System.arraycopy(w_values,0,w_ext,1,w_values.length);
		w = Utility.ArrayDivision(Utility.toObjectArray(w_ext), Utility.ArrayNorm2(Utility.toObjectArray(w_ext)));
			
		return P;
	}

	protected IloCplex initializeCplex(int Nv, int No, Vector<LinkedHashMap<Short, Double>> a_in,   Vector<Double> b_in,
			Vector<LinkedHashMap<Short, Double>> aeq_in,   Vector<Double> beq_in,double[] lb, double[] ub) throws IloException {
		this.cplex = new IloCplex();
		cplex.setWarning(null);
		cplex.setOut(null);
		if(Nv< 2000)
		{
			cplex.setParam(IloCplex.DoubleParam.DetTiLim, 100);
			cplex.setParam(IloCplex.DoubleParam.TiLim, 100);
		}
		else
		{
			cplex.setParam(IloCplex.DoubleParam.WorkMem ,2000.0);
			cplex.setParam(IloCplex.DoubleParam.DetTiLim, 10000);
			//cplex.setParam(IloCplex.DoubleParam.TiLim, 2000);
		}

		if(xVar ==null)
		{
 			String[] names = new String[Nv+No+1];
			for(int i=0;i<names.length;i++)
			{
				names[i]= "x_"+i;
			}
			xVar =  cplex.numVarArray(Nv+No+1, lb,ub, names); 
		}
		 
		// add the normal constraints, A * X <= B
		List<IloRange> constantConsts = new ArrayList<IloRange>();
		List<IloNumExpr> inEqualconsts = new ArrayList<IloNumExpr>();
		
		  // add the new inequality constraints	 
		for (int j =  0 ; j <a_in.size()  ; j++) {
			LinkedHashMap<Short, Double> array =  a_in.get(j);
			List<IloNumExpr> inEqual = new ArrayList<IloNumExpr>();
			for (Short key: array.keySet()) {
					IloNumExpr itme = cplex.prod(1.0*array.get(key) ,xVar[key]);
					inEqual.add(itme);
			}
			IloNumExpr itemsum = cplex.sum((IloNumExpr[]) inEqual.toArray(new IloNumExpr[0]));
			inEqualconsts.add(itemsum);
			IloRange extra1 = cplex.addLe(itemsum, (double) b_in.get(j));
			constantConsts.add(extra1);
		}
	

		// add the normal constraints  A_eq * X = B_eq
		List<IloNumExpr> equalconsts = new ArrayList<IloNumExpr>();
		// add the new equality constraints	 
		for (int j =  0 ; j <aeq_in.size()  ; j++) {
			LinkedHashMap<Short, Double> array =  aeq_in.get(j);
			List<IloNumExpr> equal = new ArrayList<IloNumExpr>();
			for (Short key: array.keySet()) {
					IloNumExpr itme = cplex.prod(1.0*array.get(key) , xVar[key]);
					equal.add(itme);
			}
			IloNumExpr itemsum = cplex.sum((IloNumExpr[]) equal.toArray(new IloNumExpr[0]));
			equalconsts.add(itemsum);
			IloRange extra2 = cplex.addEq(itemsum, (double) beq_in.get(j));
			constantConsts.add(extra2);
		}
	
        assert  constantConsts.size() == inEqualconsts.size()+equalconsts.size() ;
        return cplex;
	}
	
	@SuppressWarnings({ "unchecked" })
	public Vector<Double[]> calculate_old() throws Exception {
		// TODO Auto-generated method stub
		/**  Matlab code :
		  No = size(y_up,1);                % Calculate number of objectives
		    ini = rand(1,No);
		    ini = ini/sum(ini);
		    X0 = ini*y_up;                    % get an initial point on utopia plane

		    p = X0;                           % put X0 into p
		 */
		
		int No = y_up.length;
		Double[] ini = Utility.generateRandomDouble(No);
		ini = Utility.ArrayDivision(ini, Utility.arraySum(ini));
		Double[] X0 = Utility.ArrayProduceMatrix(ini, Utility.MatrixTranspose(y_up));
		Double[] p = X0;
		P.add(p);
		
		/**  Matlab code :
		 %% Find the vectors
		    V = zeros(1,No);
		    for i=1:(No-1)
		        Vi = y_up(No,:)-y_up(i,:);
		        if norm(Vi) ~=0
		            Vi = Vi/norm(Vi);
		        end;
		        V = [V;Vi];
		    end
		    V = V(2:No,:);	
		*/
		Double[][] V = new Double[No-1][No];
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
	
		/**  Matlab code :
		 %% Finding other N-1 points 
		    for i=1:(N-1)
		        mult = rand(1, No-1);           % uniformly distributed multiplier
		        mult = mult/norm(mult);
		        D = mult*V;                     % a random direction on utopia plane
		        
		        %% Finding the limit of lambda
		        lambda_l=0;
		        lambda_u=0;
		        load mo_sql_machine_web
		        %load original_problem
		        % [X,Y,lambda]
		        Nv = size(A,2);
		        A = [A,zeros(size(A,1),No+1)];
		        Aeq = [Aeq,zeros(size(Aeq,1),No+1)];
		        Aeq1 = [f,-eye(No),zeros(No,1)];
		        Aeq2 = [zeros(No,Nv),eye(No),-D'];
		        Aeq = [Aeq;Aeq1;Aeq2];
		        beq = [beq;zeros(No,1);X0'];
		        ff = [zeros(1,Nv+No),1];
		        intcon = [];
		        lb = [zeros(1,Nv),-Inf*ones(1,No),-Inf];
		        ub = [ones(1,Nv),Inf*ones(1,No),Inf];
		        [X_l,FVAL_l,exitflag_l] = intlinprog(ff,intcon,A,b,Aeq,beq,lb,ub);
		        [X_u,FVAL_u,exitflag_u] = intlinprog(-ff,intcon,A,b,Aeq,beq,lb,ub);
		        if exitflag_l==1
		            lambda_l = FVAL_l;
		        end
		        if exitflag_u==1
		            lambda_u = -FVAL_u;
		        end
		        
		        
		        %% randomly choose a new point
		        lambda = unifrnd(lambda_l,lambda_u);        
		        X0 = X0+lambda*D;
		        p = [p;X0];
		    end
		    */
		// Finding other N-1 points 
		for(int i =0; i< n-1;i++)
		{
			this.extra_A = (Vector<LinkedHashMap<Short, Double>>) ori_A.clone();
			//this.extra_B = (Vector<Double>) ori_B.clone();
			this.extra_Aeq =  (Vector<LinkedHashMap<Short, Double>>) ori_Aeq.clone();
			this.extra_Beq=  (Vector<Double>) ori_Beq.clone();
			Double[] mult = Utility.randDistributedArray(1,No-1)[0];
			mult = Utility.ArrayDivision(mult,Utility.ArrayNorm2(mult) );
			Double[] D = Utility.ArrayProduceMatrix(mult,  Utility.MatrixTranspose(V));
			//Finding the limit of lambda
			 double lambda_l=0;
			 double lambda_u=0;
			 int Nv = this.varNo;
		     //expendSparseMatrixWithZeroColumns(this.extra_A, Nv, No+1);   
		     //expendSparseMatrixWithZeroColumns(this.extra_Aeq, Nv, No+1); 
		     //modify extra_Aeq
		     
		     Vector<LinkedHashMap<Short, Double>> Aeq1 = convertMatrixs2SparseMat(f, Utility.negMatrix(Utility.eyeMatrix(No)), Utility.zeroMatrix(No,1));
		    //		 [f,-eye(No),zeros(No,1)];
		     Vector<LinkedHashMap<Short, Double>> Aeq2 = convertMatrixs2SparseMat(Utility.zeroMatrix(No,Nv),Utility.eyeMatrix(No), Utility.negMatrix(Utility.MatrixTranspose(Utility.twoDemensionize(D))));
		     //   Aeq2 = [zeros(No,Nv),eye(No),-D'];
		     this.extra_Aeq.addAll(Aeq1);
		     this.extra_Aeq.addAll(Aeq2);
		     //   Aeq = [Aeq;Aeq1;Aeq2];
		     
		     Double[] Beq1 = Utility.zeros(1, No);
		     Double[] Beq2 = X0;
		     this.extra_Beq.addAll(Arrays.asList(Beq1));
		     this.extra_Beq.addAll(Arrays.asList(Beq2));
		     //beq = [beq;zeros(No,1);X0'];
		     Double[] ff = Utility.zeros(1, Nv+No+1);
		     ff[Nv+No] = 1.0;
		     //ff = [zeros(1,Nv+No),1];
		     double[] lb = ArrayUtils.toPrimitive(Utility.zeros(1, Nv+No+1));
		     for(int k = Nv;k<lb.length; k++)
		     {
		    	 lb[k] = Double.NEGATIVE_INFINITY;
		    	 //lb = [zeros(1,Nv),-Inf*ones(1,No),-Inf];
		     }
		     double[] ub = ArrayUtils.toPrimitive(Utility.ones(1, Nv+No+1));
		     for(int k = Nv;k<ub.length; k++)
		     {
		    	 ub[k] = Double.POSITIVE_INFINITY;
				 // ub = [ones(1,Nv),Inf*ones(1,No),Inf];
		     }
		     
		     long startTime=System.currentTimeMillis();   //é‘¾å³°å½‡é”Ÿï¿½?æ¿®å¬«æ¤‚é”Ÿï¿½?  
		     CplexResult positiveRst = NCGOP.mixintlinprog (null, null, ff,this.extra_A,this.ori_B,this.extra_Aeq,this.extra_Beq, lb,ub);
		     long endTime=System.currentTimeMillis(); //é‘¾å³°å½‡ç¼�æ’´æ½«é�ƒå •æ£¿  
		     long time1= (endTime-startTime)/1000;
		     startTime=System.currentTimeMillis();   //é‘¾å³°å½‡é”Ÿï¿½?æ¿®å¬«æ¤‚é”Ÿï¿½?  
		     CplexResult negativeRst = NCGOP.mixintlinprog (null, null, Utility.negArray(ff),this.extra_A,this.ori_B,this.extra_Aeq,this.extra_Beq, lb,ub);	
		     endTime=System.currentTimeMillis(); //é‘¾å³°å½‡ç¼�æ’´æ½«é�ƒå •æ£¿
		     long time2= (endTime-startTime)/1000;
		     
		     System.out.println("for p_ié”Ÿï¿½? "+ time1+"s and "+time2+"s" );  
		     
		     if(positiveRst.getExitflag())
		     {
		    	 lambda_l= positiveRst.getFVAL();
       	  	 }
		     if(negativeRst.getExitflag())
		     {
		    	 lambda_u = negativeRst.getFVAL() *-1.0;
		     }
		     
		     double lambda = Utility.unifrnd(lambda_l,lambda_u);
		     X0 =  Utility.ArraySum(X0,Utility.ArrayMultiply(D,lambda));
		     P.add(X0);
		}
		return P;
	}

	public void setParas(Vector<LinkedHashMap<Short, Double>> ori_A2, Vector<Double> ori_B2,
			Vector<LinkedHashMap<Short, Double>> ori_Aeq2, Vector<Double> ori_Beq2,
			Double[][] f_in) {
		this.ori_A= ori_A2;
		this.ori_B = ori_B2;
		this.ori_Aeq =ori_Aeq2;
		this.ori_Beq = ori_Beq2;
		//this.extra_A = new Vector<LinkedHashMap<Short, Double>>();
		//this.extra_B = new Vector<Double>();
		this.f = f_in;
		//this.extra_Aeq = new Vector<LinkedHashMap<Short, Double>>();
		//this.extra_Beq= new Vector<Double>();
		try {
			int No = y_up.length;
			int Nv = this.varNo;
			// generate constant upper bounds and lower bounds
			this.lb = ArrayUtils.toPrimitive(Utility.zeros(1, Nv + No + 1));
			for (int k = Nv; k < lb.length; k++) {
				lb[k] = Double.NEGATIVE_INFINITY;
				// lb = [zeros(1,Nv),-Inf*ones(1,No),-Inf];
			}
			this.ub = ArrayUtils.toPrimitive(Utility.ones(1, Nv + No + 1));
			for (int k = Nv; k < ub.length; k++) {
				ub[k] = Double.POSITIVE_INFINITY;
				// ub = [ones(1,Nv),Inf*ones(1,No),Inf];
			}

			cplex = initializeCplex(Nv, No, this.ori_A, this.ori_B, this.ori_Aeq, this.ori_Beq, lb, ub);
			this.extra_A = new Vector<LinkedHashMap<Short, Double>>();
			this.extra_B = new Vector<Double>();
			this.extra_Aeq = new Vector<LinkedHashMap<Short, Double>>();
			this.extra_Beq = new Vector<Double>();
			Vector<LinkedHashMap<Short, Double>> Aeq1 = this
					.convertMatrixs2SparseMat(this.f,
							Utility.negMatrix(Utility.eyeMatrix(No)),
							Utility.zeroMatrix(No, 1));
			this.extra_Aeq.addAll(Aeq1);
			Double[] Beq1 = Utility.zeros(1, No);
			this.extra_Beq.addAll(Arrays.asList(Beq1));
			cplex = initializeCplexWihtMoreCons(this.extra_A,this.extra_B,this.extra_Aeq, this.extra_Beq);
		}
		catch(Exception e)
		{
			
		}
	}

	protected IloCplex initializeCplexWihtMoreCons(Vector<LinkedHashMap<Short, Double>> extra_A_in, Vector<Double> extra_B_in,
			Vector<LinkedHashMap<Short, Double>> extra_Aeq_in, Vector<Double> extra_Beq_in) throws IloException {
		// add the normal constraints, A * X <= B
		List<IloRange> constantConsts = new ArrayList<IloRange>();
		List<IloNumExpr> inEqualconsts = new ArrayList<IloNumExpr>();

		// add the new inequality constraints
		for (int j = 0; j < extra_A_in.size(); j++) {
			LinkedHashMap<Short, Double> array = extra_A_in.get(j);
			List<IloNumExpr> inEqual = new ArrayList<IloNumExpr>();
			for (Short key : array.keySet()) {
				IloNumExpr itme = cplex.prod(1.0 * array.get(key), xVar[key]);
				inEqual.add(itme);
			}
			IloNumExpr itemsum = cplex.sum((IloNumExpr[]) inEqual.toArray(new IloNumExpr[0]));
			inEqualconsts.add(itemsum);
			IloRange extra1 = cplex.addLe(itemsum, (double) extra_B_in.get(j));
			constantConsts.add(extra1);
		}

		// add the normal constraints A_eq * X = B_eq
		List<IloNumExpr> equalconsts = new ArrayList<IloNumExpr>();
		// add the new equality constraints
		for (int j = 0; j < extra_Aeq_in.size(); j++) {
			LinkedHashMap<Short, Double> array = extra_Aeq_in.get(j);
			List<IloNumExpr> equal = new ArrayList<IloNumExpr>();
			for (Short key : array.keySet()) {
				IloNumExpr itme = cplex.prod(1.0 * array.get(key), xVar[key]);
				equal.add(itme);
			}
			IloNumExpr itemsum = cplex.sum((IloNumExpr[]) equal.toArray(new IloNumExpr[0]));
			equalconsts.add(itemsum);
			IloRange extra2 = cplex.addEq(itemsum, (double) extra_Beq_in.get(j));
			constantConsts.add(extra2);
		}

		assert constantConsts.size() == inEqualconsts.size() + equalconsts.size();
		return this.cplex;
	}

	public static void expendSparseMatrixWithZeroColumns(Vector<LinkedHashMap<Short, Double>> ori_A, int startPos, int length) {
		 for(LinkedHashMap<Short, Double> map: ori_A)
		 {
			 for(int i=0; i< length; i++)
			 {
				 map.put((short)(startPos+i), 0.0);
			 }
		 }
	}
	

	public static Vector<LinkedHashMap<Short, Double>> convertMatrixs2SparseMat(Double[][] f2, Double[][] negMatrix,
			Double[][] zeroMatrix) {
		assert f2.length == negMatrix.length;
		assert negMatrix.length == zeroMatrix.length;
		
		Vector<LinkedHashMap<Short, Double>> sparseMapVector = new Vector<LinkedHashMap<Short, Double>>();
		List<Double[][]> matrixList = new ArrayList<Double[][]>();
		matrixList.add(f2);
		matrixList.add(negMatrix);
		matrixList.add(zeroMatrix);
		int totalLength = 0;
		for(int i=0; i< f2.length;i++)
		{
			totalLength =0;
			//for every raw build a LinkedHashMap<Short, Double>
			LinkedHashMap<Short, Double> sparseMap = new LinkedHashMap<Short, Double>();
			for(Double[][]  matrix: matrixList)
			{
				Double[] row = matrix[i];
				for(int j=0; j< row.length; j++)
				{
					double value = row[j];
					if(value!=0.0)  sparseMap.put((short) (j+totalLength), value);
				}
				totalLength+= row.length;
			}
			sparseMapVector.add(sparseMap);
		}
		return sparseMapVector;
	}
}
