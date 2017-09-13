/**
 * 
 */
package ncgop.cplex.tsl.ntu.sg;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.lang3.ArrayUtils;

import ProductLine.FeatureModel.FeatureModel;
import ProductLine.FeatureModel.LogicFeatureModel;
import cplex.tsl.ntu.sg.CplexResult;
import cplex.tsl.ntu.sg.CplexResultComparator;
import cplex.tsl.ntu.sg.CplexSolution;
import cplex.tsl.ntu.sg.MyIloCplex;
import cplex.tsl.ntu.sg.Utility;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

/**
 * @author yinxing
 *
 */
public class NCGOP_CutTry {
	public static int N = 1000; 
	public static int EXE_TIME = 0; 
	
	// Parameters Settings
	//protected MyIloCplex cplexOri;
	protected IloCplex cplexNum;
	protected MyIloCplex cplexInt;
	protected IloNumVar[] xNumVar;
	protected IloIntVar[] xIntVar;
	
	protected int objNo;   // Number of objective functions
	protected int varNv;  // Number of variables
	protected Double[][] F;
	protected Double[][] V;
	
	private List<LinkedHashMap<Short, Double>> sparseEquations;
	private List<LinkedHashMap<Short, Double>> sparseAtMostInequations;
	
	private Vector<LinkedHashMap<Short, Double>> ori_A;
	private Vector<Double> ori_B;
	
	private Vector<LinkedHashMap<Short, Double>> ori_Aeq;
	private Vector<Double> ori_Beq;
	
	private Vector<LinkedHashMap<Short, Double>> extra_A;
	private Vector<Double> extra_B;
	private Vector<LinkedHashMap<Short, Double>> extra_Aeq;
	private Vector<Double> extra_Beq;
	
	// InequationMap inequationMap;
	protected int omittedBits;
	
	public Double[][] getF() {
		return F;
	}
    
	//new field added to support the cut&try approach
	private Double[] x;
	private Double[] y;
	private Double[] yy;
	private Double [] yyy;
	private Vector<Double[]> nonDominantSols;
	
	protected Set<Boolean[]> E_out;

	private UtopiaPlane utopiaPlane;
	private AlphaLimit alphaLimit;
	private SolRep_CutTry pGenerator;

	public UtopiaPlane getUtopiaPlane() {
		return utopiaPlane;
	}


	public AlphaLimit getAlphaLimit() {
		return alphaLimit;
	}

	public SolRep_CutTry getpGenerator() {
		return pGenerator;
	}

	public Set<Boolean[]> getE_out() {
		return E_out;
	}

	public Vector<LinkedHashMap<Short, Double>> getOri_A() {
		return ori_A;
	}


	public Vector<Double> getOri_B() {
		return ori_B;
	}


	public Vector<LinkedHashMap<Short, Double>> getOri_Aeq() {
		return ori_Aeq;
	}


	public Vector<Double> getOri_Beq() {
		return ori_Beq;
	}


	public Vector<LinkedHashMap<Short, Double>> getExtra_A() {
		return extra_A;
	}
 

	public Vector<Double> getExtra_B() {
		return extra_B;
	}
	
	public NCGOP_CutTry(MyIloCplex cplex, int objNo, int varNv) throws Exception {
		this.objNo = objNo;
		this.varNv= varNv;
		FeatureModel fm = cplex.getNewFeatureModel() ;
		if(fm instanceof LogicFeatureModel)
		{
			LogicFeatureModel lfm = (LogicFeatureModel)fm;
			omittedBits = lfm.getMandatoryNames().size()+ lfm.getMustNotInNames().size();
		}
		else
		{
			omittedBits  = 0; 
			throw new Exception("input model not LogicFeatureModel");
		}
		E_out = new LinkedHashSet<Boolean[]>();
		this.cplexInt = cplex;
		this.nonDominantSols = new Vector<Double[]>();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {

	}


	public void setParas(Double[][] f2, List<LinkedHashMap<Short, Double>> sparseEquations2,
			List<LinkedHashMap<Short, Double>> sparseAtMostInequations2) {
		this.F = f2;
	    this.sparseEquations =sparseEquations2;
	    this.sparseAtMostInequations = sparseAtMostInequations2;
		
		this.extra_A = new Vector<LinkedHashMap<Short, Double>>();
		//this.extra_A.addAll(Arrays.asList(this.ori_A));
	    this.extra_B = new Vector<Double>( ); 
	    //this.extra_B.addAll(Arrays.asList(this.ori_B));
	    
	    //add contraints of  sparseEquations, sparseAtMostInequations to the cplex
	    try {
	    	builtEqualAndInequalMaps();
			initializeCplexInt();
			initializeCplexNum();
		} catch (IloException e) {
			e.printStackTrace();
		}
		
	}

	@SuppressWarnings("unchecked")
	private void builtEqualAndInequalMaps() {
		// TODO Auto-generated method stub
		this.ori_A = new Vector<LinkedHashMap<Short, Double>>();
		this.ori_B = new Vector<Double>();
		
		this.ori_Aeq = new Vector<LinkedHashMap<Short, Double>>();
		this.ori_Beq = new Vector<Double>();
		
		for(int i = 0; i<this.sparseAtMostInequations.size();i++ )
		{
			LinkedHashMap<Short, Double> inequal = this.sparseAtMostInequations.get(i);
			LinkedHashMap<Short, Double> newinequal =  (LinkedHashMap<Short, Double>) inequal.clone();
			Double right = newinequal.remove((short)this.varNv);
			this.ori_A.addElement(newinequal);
			this.ori_B.addElement(right);
		}
		
		for(int i = 0; i<this.sparseEquations.size();i++ )
		{
			LinkedHashMap<Short, Double> equal = this.sparseEquations.get(i);
			LinkedHashMap<Short, Double> newequal =  (LinkedHashMap<Short, Double>) equal.clone();
			Double right = newequal.remove((short)this.varNv);
			this.ori_Aeq.addElement(newequal);
			this.ori_Beq.addElement(right);
		}
		
		assert 	this.ori_A.size() == this.ori_B.size();
		assert 	this.ori_Aeq.size() == this.ori_Beq.size();
	}


	protected void initializeCplexInt() throws IloException {
		//this.cplexInt= new IloCplex();
		cplexInt.setOut(null);
		cplexInt.setWarning(null);
 
		List<String> varNames = new ArrayList<String>();
		for(int i = 1; i<= this.varNv; i++)
		{
			varNames.add("X_"+i);
			
		}
		IloIntVar[] x = cplexInt.intVarArray(varNv, 0, 1, varNames.toArray(new String[0]));
		this.xIntVar = x;
		// add the normal constraints, A * X <= B
		List<IloNumExpr> inEqualconsts = new ArrayList<IloNumExpr>();
		List<IloRange> constantConsts = new ArrayList<IloRange>();
		for (LinkedHashMap<Short, Double> array: this.sparseAtMostInequations) {
			//Byte[] array = ori_A[j];
			List<IloNumExpr> inEqual = new ArrayList<IloNumExpr>();
			for (Short key: array.keySet()) {
				    if(key == varNv) continue;
					IloNumExpr itme = cplexInt.prod(1.0*array.get(key) , x[key]);
					inEqual.add(itme);
			}
			IloNumExpr itemsum = cplexInt.sum((IloNumExpr[]) inEqual.toArray(new IloNumExpr[0]));

			inEqualconsts.add(itemsum);
			IloRange inequalConst = cplexInt.addLe(itemsum, array.get((short)varNv));
			constantConsts.add(inequalConst);
		}

		// add the normal constraints  A_eq * X = B_eq
		List<IloNumExpr> equalconsts = new ArrayList<IloNumExpr>();
		for (LinkedHashMap<Short, Double> array: this.sparseEquations) {
			//Byte[] array = this.A_eq[j];
			List<IloNumExpr> equal = new ArrayList<IloNumExpr>();
			for (Short key: array.keySet()) {
				    if(key == varNv) continue;
					IloNumExpr itme = cplexInt.prod(1.0 * array.get(key), x[key]);
					equal.add(itme);
			}
			IloNumExpr itemsum = cplexInt.sum((IloNumExpr[]) equal.toArray(new IloNumExpr[0]));
			// IloNumExpr inEqualconst =
			// cplex.sum(itemsum,array[array.length-1]);
			equalconsts.add(itemsum);
			IloRange equalConst =  cplexInt.addEq(itemsum,array.get((short)varNv));
			constantConsts.add(equalConst);
		}
	
        assert  constantConsts.size() == this.sparseEquations.size()+this.sparseAtMostInequations.size() ;
	}

	protected void initializeCplexNum() throws IloException {
		this.cplexNum= new IloCplex();
		cplexNum.setOut(null);
		cplexNum.setWarning(null);
 
		List<String> varNames = new ArrayList<String>();
		for(int i = 1; i<= this.varNv+1; i++)
		{
			varNames.add("x_"+i);
			
		}
		double[] lb = Utility.toPrimateArray(Utility.zeros(1, varNv+1));
		lb[varNv]=  Double.NEGATIVE_INFINITY;
		double[] ub = Utility.toPrimateArray(Utility.ones(1, varNv+1));
		ub[varNv]= 0.0;
		IloNumVar[] x = cplexNum.numVarArray(varNv+1, lb, ub, varNames.toArray(new String[0]));
		xNumVar = x;
		
		// add the normal constraints, A * X <= B
		List<IloNumExpr> inEqualconsts = new ArrayList<IloNumExpr>();
		List<IloRange> constantConsts = new ArrayList<IloRange>();
		for (LinkedHashMap<Short, Double> array: this.sparseAtMostInequations) {
			//Byte[] array = ori_A[j];
			List<IloNumExpr> inEqual = new ArrayList<IloNumExpr>();
			for (Short key: array.keySet()) {
				    if(key == varNv) continue;
					IloNumExpr itme = cplexNum.prod(1.0*array.get(key) , x[key]);
					inEqual.add(itme);
			}
			IloNumExpr itemsum = cplexNum.sum((IloNumExpr[]) inEqual.toArray(new IloNumExpr[0]));

			inEqualconsts.add(itemsum);
			IloRange inequalConst = cplexNum.addLe(itemsum, array.get((short)varNv));
			constantConsts.add(inequalConst);
		}

		// add the normal constraints  A_eq * X = B_eq
		List<IloNumExpr> equalconsts = new ArrayList<IloNumExpr>();
		for (LinkedHashMap<Short, Double> array: this.sparseEquations) {
			//Byte[] array = this.A_eq[j];
			List<IloNumExpr> equal = new ArrayList<IloNumExpr>();
			for (Short key: array.keySet()) {
				    if(key == varNv) continue;
					IloNumExpr itme = cplexNum.prod(1.0 * array.get(key), x[key]);
					equal.add(itme);
			}
			IloNumExpr itemsum = cplexNum.sum((IloNumExpr[]) equal.toArray(new IloNumExpr[0]));
			// IloNumExpr inEqualconst =
			// cplex.sum(itemsum,array[array.length-1]);
			equalconsts.add(itemsum);
			IloRange equalConst =  cplexNum.addEq(itemsum,array.get((short)varNv));
			constantConsts.add(equalConst);
		}
	
        assert  constantConsts.size() == this.sparseEquations.size()+this.sparseAtMostInequations.size() ;
	}
	
	public void execute(Vector<LinkedHashMap<Short, Double>> a_in, Vector<Double> b_in, Double[][] f_in, Double[][] f_out, int k,  Set<Boolean[]> e_in,Map<String, CplexSolution> sols ) throws Exception
	{
		//for utopiaPlane calculation, to assure the completeness of resolving, we cannot set timeout for this.cplex. 
		//then for NCPDG, we can give the timeout for intlinprog timeout for each execution 
		// for non-Linux projects
		if(this.varNv< 2000)
		{
			cplexInt.setParam(IloCplex.DoubleParam.DetTiLim, 100);
			cplexInt.setParam(IloCplex.DoubleParam.TiLim, 100);
			cplexNum.setParam(IloCplex.DoubleParam.DetTiLim, 100);
			cplexNum.setParam(IloCplex.DoubleParam.TiLim, 100);
		}
		else
		{
			cplexInt.setParam(IloCplex.DoubleParam.WorkMem ,2000.0);
			cplexInt.setParam(IloCplex.DoubleParam.DetTiLim, 10000);
			cplexNum.setParam(IloCplex.DoubleParam.WorkMem ,2000.0);
			cplexNum.setParam(IloCplex.DoubleParam.DetTiLim, 10000);
//			cplex.setParam(IloCplex.DoubleParam.DetTiLim, 5000);
//			cplex.setParam(IloCplex.DoubleParam.TiLim, 5000);
		}
		
		this.utopiaPlane = new UtopiaPlane(this.cplexInt,this.xIntVar,a_in,b_in,f_in);
		this.utopiaPlane.calculate();
		System.out.println("utopiaPlane done.");
	 
		this.pGenerator = new SolRep_CutTry(this.utopiaPlane.getY_up(), this.varNv,N);
		this.pGenerator.setParas(this.ori_A,this.ori_B,this.ori_Aeq,this.ori_Beq,f_in);
		Vector<Double[]> points = this.pGenerator.calculateWandV();
		System.out.println("w, V calculation done.");
	    //this.V = calculateV(this.utopiaPlane.getY_up(),this.objNo);
		int n_sol=0;
		Double[] p_k = pGenerator.P.get(0);
		int counter =0; 
		while(n_sol<N+1){
			if(n_sol!=0)
			{
				p_k= HitAndRunGenerator.calculateReferPoint(p_k,pGenerator);
				//p_k= Utility.dummyPkForTesting_0911a2();
				pGenerator.P.add(p_k);
			}
			
			System.out.println("using p_k: "+counter++);
			CplexResult result = calculate(this.F,this.ori_A, this.ori_B, this.ori_Aeq, this.ori_Beq,this.utopiaPlane.getY_up(),p_k,this.utopiaPlane.getY_ub(), this.utopiaPlane.getY_lb(),this.pGenerator.w,ConstantMatrix.v, sols);
			if(result!=null && checkNonDominance(result, sols))
			{
				addSolutionXToMap(result, sols);
			}
			n_sol= sols.size();
		}
		
	    //int counter =0; 
		//for(Double[] p_k: points)
		//{
		//	System.out.println("using p_k: "+counter++);
		//	calculate(this.F,this.ori_A, this.ori_B, this.ori_Aeq, this.ori_Beq,this.utopiaPlane.getY_up(),p_k,this.utopiaPlane.getY_ub(), this.utopiaPlane.getY_lb(),sols);
		//}
		System.out.println("Find solution num: "+sols.size());
	}

	private boolean checkNonDominance(CplexResult result, Map<String, CplexSolution> sols) {
		// TODO Auto-generated method stub
		String key =extractSolution (result); 
		if(sols.containsKey(key))
			return false;
		else
		{		
			boolean nonDominant = checkNonDominance(key, this.nonDominantSols);
			return nonDominant;
		}
	}

	public static boolean checkNonDominance(String key, Vector<Double[]> nonDominantSols2) {
		String[] values = key.split("_");
		Double[] doubleValues = new Double[values.length];
		for(int i=0; i < values.length; i++)
		{
			doubleValues[i] = Double.parseDouble(values[i]);
		}
		
		for(Double[] currentSol: nonDominantSols2)
		{
			if(Utility.dominateSol(currentSol,doubleValues))
			{
				return false;
			}
		}
		nonDominantSols2.add(doubleValues);
		return true;
	}

	public CplexResult calculate(Double[][] f,
			Vector<LinkedHashMap<Short, Double>> ori_A2, Vector<Double> ori_B2,
			Vector<LinkedHashMap<Short, Double>> ori_Aeq2,
			Vector<Double> ori_Beq2, Double[][] y_up, Double[] p_k,
			Double[] y_ub, Double[] y_lb, Double[] w, Double[][] v, Map<String, CplexSolution> sols) throws Exception {
 
		CplexResult returnedResult = null;
		int No = this.objNo;
		int Nv = this.varNv;
	 
		/**  Matlab code :
		%add constraint: f*x = p_k + lambda*w
			AAeq = [Aeq,zeros(size(Aeq,1),1)];
			AAeq = [AAeq;[f,-w]];
			bbeq = [beq;p_k];
			AA = [A,zeros(size(A,1),1)];
		*/
	    //[Aeq,zeros(size(Aeq,1),1)] is to add one column of 0 for the existing Aeq
		//get the [f,-w]
	    Double[] neg_w = Utility.negArray(w);
	    //get the extra Aeq
	    this.extra_Aeq = Utility.denseTwoMatrix2SparseMatrix(f, Utility.twoDemensionizeAndTranspose(neg_w));
	    Vector<LinkedHashMap<Short, Double>> tempAeq=  this.extra_Aeq;
	    //get the extra Beq
	    this.extra_Beq = Utility.denseArray2SparseArray(p_k);
	    Vector<Double> tempBeq = this.extra_Beq;
		
	    //while no extra for the inequalities A and B
		
	    /**  Matlab code :
	    lb = [zeros(1,Nv),-Inf];                                      % Lower and upper bound of variables
	    ub = [ones(1,Nv),0] ;
	    intcon = [];

	    ff = [zeros(1,Nv),1]; 

	    [X,FVAL,Exitflag] = intlinprog (ff,intcon,AA,b,AAeq,bbeq,lb,ub);
	    */
	    double[] lb = Utility.toPrimateArray(Utility.zeros(1, Nv  + 1));
		lb[Nv] = Double.NEGATIVE_INFINITY;
		 
		double[] ub = Utility.toPrimateArray(Utility.ones(1, Nv + 1));
		ub[Nv] = 0.0;

		Double[] lb_Obj =(Utility.zeros(1, Nv  + 1));
		lb_Obj[Nv] = Double.NEGATIVE_INFINITY;
			 
		Double[] ub_Obj =(Utility.ones(1, Nv + 1));
		ub_Obj[Nv] = 0.0;
			
	    Double[] ff = Utility.zeros(1, Nv  + 1);
	    ff[Nv] = (double) 1; 
	    
		CplexResult result1 = NCGOP_CutTry.mixintlinprog(cplexNum,
				this.xNumVar, ff, null, null, tempAeq, tempBeq,
				lb, ub);
		
		if(result1.getExitflag() == true)
		{
			/**  Matlab code :
			yy = (f*X(1:Nv,:))';
			%% solving the 2nd problem, vertical constraint
			AA = [A;f(1:No-1,:)];                                        % add constraints
			bb = [b;yy(:,1:No-1)'];
			intcon = 1:Nv;
			% calculate new objective function with CWMOIP
			ff = f(No,:);
			*/
			Double[] X= new Double[Nv];
			System.arraycopy(result1.getXvar(), 0, X, 0, X.length); //get the values for the Nv deciding variables
			this.yy = Utility.ArrayProduceMatrix(X, f); 
			//f(1:No-1,:) is to get the first No-1 rows of objective expression 
		    Double[][] AA = Utility.getFirstItem(f, No-1);
		    Vector<LinkedHashMap<Short, Double>> AAtemp = Utility.denseMatrix2SparseMatrix(AA); // the extra_A is the f(1:No-1,:)
		    // the extra_b is the yy(:,1:No-1)'
		    Double[] bb = Utility.getFirstItem(yy, No-1);
		    Vector<Double> bbtemp = Utility.denseArray2SparseArray(bb);
			// intcon = 1:Nv --- all the deciding variables are int, int programming
			// ff = f(No,:). the last row (objective) of ff
		    ff = f[No-1];
		    /**  Matlab code :
		    for i=1:(No-1)
		            i = No-i;
		            weight = 1/(y_ub(i,1)-y_lb(i,1)+1);
		            ff = ff + weight*f(i,:);
		    end
		    */
		    //getting the constraint-weighted ff
		    for(int i=No-2; i>= 0; i--){
		    	int pos=i;
		    	double weight = 1.0/(y_ub[pos]-y_lb[pos]+1);
		    	int targetPos =No-2- i;
		    	Double[] ffDelta = Utility.ArrayMultiply(f[pos],weight);
		    	ff = Utility.ArraySum(ff, ffDelta);
		    }
		    /**  Matlab code :
		     *  [X,FVAL,Exitflag2] = intlinprog (ff,intcon,AA,bb,Aeq,beq,lb,ub);
		     */
		    CplexResult result2 = NCGOP_CutTry.intlinprog(cplexInt,
					this.xIntVar, ff, AAtemp, bbtemp, null, null,
					lb_Obj, ub_Obj);
		    
		    if(result2.getExitflag()!= true) //% solving until find a solution
		    {
		    	 /**  Matlab code :
		    		        bb = [b;p_k(1:No-1,1)];
		    		        [X,FVAL,Exitflag3] = intlinprog (ff,intcon,AA,bb,Aeq,beq,lb,ub);
		    		        if Exitflag3 ==1
		    		            p = p_k;
		    	  */
		    	bb = Utility.getFirstItem(p_k, No-1);
		    	bbtemp = Utility.denseArray2SparseArray(bb);
		    	CplexResult result3 = NCGOP_CutTry.intlinprog(cplexInt,
						this.xIntVar, ff, AAtemp, bbtemp, null, null,
						lb_Obj, ub_Obj);
		    	
		    	if(result3.getExitflag()==true)
		    	{
		    		Double[] p = p_k;
		    		/**  Matlab code :
		    		  while Exitflag3 ==1
	    		                X = round(X);
	    		                y = (f*X)';
	    		                x = X';
	    		                p = 1/2*(p+yy');
	    		                bb = [b;p(1:No-1,1)];
	    		                [X,FVAL,Exitflag3] = intlinprog (ff,intcon,AA,bb,Aeq,beq,lb,ub);
	    		            end
	    		     */
		    		while(result3.getExitflag()==true)
		    		{
		    			X= result3.getXvar();
		    			y =  evaluateOnMultiObjs(result3,f);
		    			returnedResult= result3;
		    			this.x= X;
		    			p = Utility.ArrayMultiply(Utility.ArraySum(p, yy),0.5);
		    			bb = Utility.getFirstItem(p, No-1);
				    	bbtemp = Utility.denseArray2SparseArray(bb);
		    			result3 = NCGOP_CutTry.intlinprog(cplexInt,
								this.xIntVar, ff, AAtemp, bbtemp, null, null,
								lb_Obj, ub_Obj);
		    		}
		    	}
		    	else{
		    		 /**  Matlab code :
		    		 else
    		            y = zeros(1,No);
    		            x = zeros(1,Nv);
    		         */
		    		y= Utility.zeros(1, No);
		    		x= Utility.zeros(1, Nv);
		    	}//end of result3
		    } 
		    else
		    {
		    	 /**  Matlab code :
		    	  X = round(X);
		          y = (f*X)';
		          x = X';
		         */
		    	//addSolutionXToMap(result2,sols);
		    	X= result2.getXvar();
		    	y= evaluateOnMultiObjs(result2,f);
		    	x= X;
		    	returnedResult=result2;
		    }//end of result2
 		}
		else {
			/**  Matlab code :
			y = zeros(1,No);
			*/
			y= Utility.zeros(1, No);
			//no results added to sols
		} //end of result1
		/**  Matlab code :
		% epsilon constraint
		AA = [A;f(1:No-1,:)];                                        % add constraints
		bb = [b;p_k(1:No-1,:)];
		intcon = 1:Nv;
		[X,FVAL,Exitflag3] = intlinprog (ff,intcon,AA,bb,Aeq,beq,lb,ub);
		if Exitflag3==1
		    yyy = (f*X)';
		else
		    yyy = zeros(1,No);
		end
		*/
		
		//calculate the non-gradual approach
//		Double[][] AA = Utility.getFirstItem(f, No-1);
//	    Vector<LinkedHashMap<Short, Double>> AAtemp = Utility.denseMatrix2SparseMatrix(AA); // the extra_A is the f(1:No-1,:)
//	    Double[] bb = Utility.getFirstItem(p_k, No-1);
//	    Vector<Double> bbtemp = Utility.denseArray2SparseArray(bb);
//	    CplexResult result4 = NCGOP_CutTry.intlinprog(cplexInt,
//				this.xIntVar, ff, AAtemp, bbtemp, null, null,
//				lb_Obj, ub_Obj);
//		if(result4.getExitflag()==true)
//		{
//			yyy= evaluateOnMultiObjs(result4,f);
//		}
//		else {
//			yyy= Utility.zeros(1, No);
//		}
		return returnedResult;
	}

	protected void addSolutionXToMap(CplexResult rslt, Map<String, CplexSolution> sols) {
		// f[0] = this.costEffiecient;
		// f[1] = this.varietyEffiecient;
		// f[2] = this.usedbeforeEffiecient;
		// f[3] = this.defectEffiecient; 
		Double[] xvar = rslt.getXvar();
		double objval = rslt.getFVAL();
		int missingFeatureSize = (int)(xvar.length- omittedBits+ Utility.ArrayProducts(xvar,F[1], xvar.length- omittedBits));
		int notUsedFeatureSize = (int) (Utility.ArrayProducts(xvar,F[2],xvar.length- omittedBits));
		double defects = Utility.ArrayProducts(xvar,F[3], xvar.length- omittedBits);
		double costs = Utility.ArrayProducts(xvar,F[0], xvar.length- omittedBits);
		
		CplexSolution sol = new CplexSolution(CplexResultComparator.formatDouble2(costs),missingFeatureSize,notUsedFeatureSize,CplexResultComparator.formatDouble2(defects),xvar);
		//CplexSolution sol = new CplexSolution(objval,missingFeatureSize,notUsedFeatureSize,defects,xvar); 
		if(!sols.containsKey(sol.getSolutionID())) 
		{
			sols.put(sol.getSolutionID(), sol);
		}
		System.out.println("BestObj: " + objval + " Cost: "+ costs);
		//System.out.println("Var: " + Arrays.asList(xvar));
		
	}
	
	protected String extractSolution(CplexResult rslt) {
		// f[0] = this.costEffiecient;
		// f[1] = this.varietyEffiecient;
		// f[2] = this.usedbeforeEffiecient;
		// f[3] = this.defectEffiecient; 
		Double[] xvar = rslt.getXvar();
		int missingFeatureSize = (int)(xvar.length- omittedBits+ Utility.ArrayProducts(xvar,F[1], xvar.length- omittedBits));
		int notUsedFeatureSize = (int) (Utility.ArrayProducts(xvar,F[2],xvar.length- omittedBits));
		double defects = Utility.ArrayProducts(xvar,F[3], xvar.length- omittedBits);
		double costs = Utility.ArrayProducts(xvar,F[0], xvar.length- omittedBits);
		
		CplexSolution sol = new CplexSolution(CplexResultComparator.formatDouble2(costs),missingFeatureSize,notUsedFeatureSize,CplexResultComparator.formatDouble2(defects),xvar);
	    return sol.getSolutionID(); 
	}
	
	protected Double[] evaluateOnMultiObjs(CplexResult rslt, Double[][] f)
	{
		Double[] xvar = rslt.getXvar();
		Double[] evals= Utility.MatrixProduceArray(f, xvar);
		return evals;
	}

	public static List<CplexResult> intlinprog(MyIloCplex cplex, IloIntVar[] xVar,  Double[] doubles,
			Vector<LinkedHashMap<Short, Double>> a_in, Vector<Double> b_in, Double[] lb, Double[] ub) throws IloException {
		// doubles as the  coefficient for target function
		IloNumExpr[] objExp = new IloNumExpr[doubles.length];
		IloNumVar[] x = xVar;
		for (int i = 0; i < x.length; i++) {
			objExp[i] = cplex.prod(doubles[i], x[i]);
		}
		IloNumExpr expr = cplex.sum(objExp);
		IloObjective temObj  = cplex.addMinimize(expr);
		
		List<IloRange> tempConsts= new ArrayList<IloRange>();
	    // add the new inequality constraints	 
		for (int j =  0 ; j <a_in.size()  ; j++) {
			LinkedHashMap<Short, Double> array =  a_in.get(j);
			List<IloNumExpr> inEqual = new ArrayList<IloNumExpr>();
			for (Short key: array.keySet()) {
					IloNumExpr itme = cplex.prod(1.0*array.get(key) , x[key]);
					inEqual.add(itme);
			}
			IloNumExpr itemsum = cplex.sum((IloNumExpr[]) inEqual.toArray(new IloNumExpr[0]));
			IloRange extra1 = cplex.addLe(itemsum, (double) b_in.get(j));
			tempConsts.add(extra1);
		}

//		boolean exitFlag = cplex.solve(); 
		/*
	       * Since the objective value is integer, setting the absolute pool gap
	       * to 0.5 ensures that a solution worse than the incumbent (which would
	       * require its objective value to be lower than the incumbent's by 1)
	       * will be rejected. To find all feasible solutions, this and the
	       * relative pool gap (SolnPoolGap) would both be kept at their default
	       * (maximum) values.
	       */
	      cplex.setParam(IloCplex.DoubleParam.SolnPoolAGap, 0.5);
	      /*
	       * Set the pool intensity to 4 to find "all" solutions.
	       */
	      cplex.setParam(IloCplex.IntParam.SolnPoolIntensity, 4);
	      /*
	       * Set the pool capacity to hold the maximum number of solutions
	       * you really want.
	       */
	      cplex.setParam(IloCplex.IntParam.SolnPoolCapacity, 5);
	      /*
	       * Set a limit on the number of solutions to find; make it very
	       * large to find all optimal solutions. The limit needs to be larger
	       * than the pool capacity to allow for suboptimal solutions found
	       * along the way.
	       */
	      cplex.setParam(IloCplex.IntParam.PopulateLim, 10);
	      /*
	       * Set the replacement strategy to turf out suboptimal solutions
	       * if space fills up; the default is FIFO, which could keep suboptimal
	       * solutions and turf out optimal ones.
	       */
	      cplex.setParam(IloCplex.IntParam.SolnPoolReplace, 1);
	      boolean exitFlag =  cplex.populate();
	      List<CplexResult> results = new ArrayList<CplexResult>();
		//CplexResult result = null;
		if (exitFlag) {
//			double objval = cplex.getObjValue(); // ��ȡ�����е����о��߱����Ľ�ֵ��
//			double[] xval = cplex.getValues(x);
//			Double[] xvar = Utility.toObjectArray(xval);
//			result = new CplexResult(objval, xvar, exitFlag);
			results=showSolutions(cplex, x ); 
        }
		else
		{
			//result = new CplexResult(Integer.MAX_VALUE, null, exitFlag);
		}
		
		//remove the temp constraints 
        for(IloRange  extraConst: tempConsts)
        {
        	cplex.delete(extraConst);
        }
    	//remove the temp objective 
        cplex.delete(temObj);
		return results;
	}
	 
	  private static final double TOL = 1e-5;
 
	 /**
	   * Display all optimal solutions found, screening out the suboptimal ones
	   * (if any).
	 * @param cplex 
	 * @param x 
	   */
	  private static List<CplexResult> showSolutions(MyIloCplex cplex, IloNumVar[] x) throws IloException {
	    // Get the number of solutions in the pool.
	    int nsol = cplex.getSolnPoolNsolns();
	    // Create a container for the indices of optimal solutions.
	    Set<Integer> opt = new LinkedHashSet<>();
	    double best = Double.NEGATIVE_INFINITY;  // best objective value found
	    /*
	     * Check which pool solutions are truly optimal; if the pool capacity
	     * exceeds the number of optimal solutions, there may be suboptimal
	     * solutions lingering in the pool.
	     */
	    for (int i = 0; i < nsol; i++) {
	      // Get the objective value of the i-th pool solution.
	      double z = cplex.getObjValue(i);
	      if (z > best + TOL) {
	        /*
	         * If this solution is better than the previous best, the previous
	         * solutions must have been suboptimal; drop them all and count this one.
	         */
	        best = z;
	        opt.clear();
	        opt.add(i);
	      } else if (z > best - TOL) {
	        /*
	         * If this solution is within rounding tolerance of optimal, count it.
	         */
	        opt.add(i);
	      }
	    }
	    System.out.println("\n\nFound " + nsol + " solutions, of which "
	                       + opt.size() + " are optimal.");
	    
	    List<CplexResult> results = new ArrayList<CplexResult>();
	    /*
	     * Now extract the optimal solutions from the pool and print them.
	     * Note: the solution pool uses zero-based indexing.
	     */
	    int h = 1;    // cumulative index of optimal solutions
	    for (int k : opt) {  // for each index of an optimal solution ...
	      System.out.println("Solution #" + (h++) + " (obj value = " + cplex.getObjValue(k) + "):");
	      double objval = cplex.getObjValue(k); // ��ȡ�����е����о��߱����Ľ�ֵ��
	      double[] xval = cplex.getValues(x,k);
	      Double[] xvar = Utility.toObjectArray(xval);
	      CplexResult result = new CplexResult(objval, xvar, true);   
	      results.add(result);
	    }
		return results;
	  }
	  
	/**
	 * 
	 * @param cplex
	 * @param xVar
	 * @param doubles
	 * @param a_in
	 * @param b_in
	 * @param aeq_in
	 * @param beq_in
	 * @param lb
	 * @param ub
	 * @return
	 * @throws IloException
	 */
	public static CplexResult intlinprog(IloCplex cplex, IloIntVar[] xVar,  Double[] doubles,
			Vector<LinkedHashMap<Short, Double>> a_in, Vector<Double> b_in, 
			Vector<LinkedHashMap<Short, Double>> aeq_in, Vector<Double> beq_in, Double[] lb, Double[] ub) throws IloException {
		IloCplex origiCplex = cplex;
		if(cplex == null)
		{
			try {
				cplex = new IloCplex();
				cplex.setWarning(null);
				cplex.setOut(null);
				cplex.setParam(IloCplex.DoubleParam.DetTiLim, 100);
				cplex.setParam(IloCplex.DoubleParam.TiLim, 100);
				// create model and solve it
				} catch (IloException e) {
				System.err.println("Concert exception caught: " + e);
				}
		}
		if(xVar ==null|| xVar.length !=doubles.length)
		{
//			IloIntVar[] intVarArray(int n,
//                    int[] min,
//                    int[] max,
//                    java.lang.String[] name)
			String[] names = new String[doubles.length];
			for(int i=1;i<=names.length;i++)
			{
				names[i-1]= "X_"+i;
			}
			xVar =  cplex.intVarArray(doubles.length, Utility.ArrayGetFloor(lb),Utility.ArrayGetCeiling(ub), names);
		}
		
		// doubles as the  coefficient for target function
		IloNumExpr[] objExp = new IloNumExpr[doubles.length];
		IloNumVar[] x = xVar;
		for (int i = 0; i < x.length; i++) {
			objExp[i] = cplex.prod(doubles[i], x[i]);
		}
		IloNumExpr expr = cplex.sum(objExp);
		IloObjective temObj  = cplex.addMinimize(expr);
		
		List<IloRange> tempConsts= new ArrayList<IloRange>();
	    // add the new inequality constraints	 
		for (int j =  0 ; j <a_in.size()  ; j++) {
			LinkedHashMap<Short, Double> array =  a_in.get(j);
			List<IloNumExpr> inEqual = new ArrayList<IloNumExpr>();
			for (Short key: array.keySet()) {
					IloNumExpr itme = cplex.prod(1.0*array.get(key) , x[key]);
					inEqual.add(itme);
			}
			IloNumExpr itemsum = cplex.sum((IloNumExpr[]) inEqual.toArray(new IloNumExpr[0]));
			IloRange extra1 = cplex.addLe(itemsum, (double) b_in.get(j));
			tempConsts.add(extra1);
		}
		// add the new equality constraints	 
		if (aeq_in != null) {
			for (int j = 0; j < aeq_in.size(); j++) {
				LinkedHashMap<Short, Double> array = aeq_in.get(j);
				List<IloNumExpr> equal = new ArrayList<IloNumExpr>();
				for (Short key : array.keySet()) {
					IloNumExpr itme = cplex.prod(1.0 * array.get(key), x[key]);
					equal.add(itme);
				}
				IloNumExpr itemsum = cplex.sum((IloNumExpr[]) equal.toArray(new IloNumExpr[0]));
				IloRange extra2 = cplex.addEq(itemsum, (double) beq_in.get(j));
				tempConsts.add(extra2);
			}
		}
		boolean exitFlag = cplex.solve(); 

		CplexResult result = null;
		if (exitFlag) {
			double objval = cplex.getObjValue(); // ��ȡ�����е����о��߱����Ľ�ֵ��
			double[] xval = cplex.getValues(x);
			Double[] xvar = Utility.toObjectArray(xval);
			result = new CplexResult(objval, xvar, exitFlag);
        }
		else
		{
			result = new CplexResult(Integer.MAX_VALUE, null, exitFlag);
		}
				
        if(origiCplex==null)
        {
        	cplex.end();
        }
        else
        {
        	//remove the temp constraints 
            for(IloRange  extraConst: tempConsts)
            {
            	cplex.delete(extraConst);
            }
        	//remove the temp objective 
            cplex.delete(temObj);
        }
		return result;
	}
	 
	/**
	 * 
	 * @param cplex
	 * @param xVar
	 * @param doubles
	 * @param a_in
	 * @param b_in
	 * @param aeq_in
	 * @param beq_in
	 * @param lb
	 * @param ub
	 * @return
	 * @throws IloException
	 */
		public static CplexResult mixintlinprog(IloCplex cplex, IloNumVar[] xVar,  Double[] doubles,
				Vector<LinkedHashMap<Short, Double>> a_in, Vector<Double> b_in, 
				Vector<LinkedHashMap<Short, Double>> aeq_in, Vector<Double> beq_in, double[] lb, double[] ub) throws IloException {
			IloCplex origiCplex = cplex;
			if(cplex == null)
			{
				try {
					cplex = new IloCplex();
					cplex.setWarning(null);
					cplex.setOut(null);
					cplex.setParam(IloCplex.DoubleParam.DetTiLim, 100);
					cplex.setParam(IloCplex.DoubleParam.TiLim, 100);
					// create model and solve it
					} catch (IloException e) {
					System.err.println("Concert exception caught: " + e);
					}
			}
			if(xVar ==null || xVar.length !=doubles.length)
			{
				String[] names = new String[doubles.length];
				for(int i=1;i<=names.length;i++)
				{
					names[i-1]= "X_"+i;
				}
				xVar =  cplex.numVarArray(doubles.length, lb,ub, names); 
				//xVar =  cplex.intVarArray(doubles.length, Utility.ArrayGetFloor(lb),Utility.ArrayGetCeiling(ub), names);
			}
			
			// doubles as the  coefficient for target function
			IloNumExpr[] objExp = new IloNumExpr[doubles.length];
			IloNumVar[] x = xVar;
			for (int i = 0; i < x.length; i++) {
				objExp[i] = cplex.prod(doubles[i], x[i]);
			}
			IloNumExpr expr = cplex.sum(objExp);
			IloObjective temObj  = cplex.addMinimize(expr);
			
			List<IloRange> tempConsts= new ArrayList<IloRange>();
		    // add the new inequality constraints
		if (a_in != null) {
			for (int j = 0; j < a_in.size(); j++) {
				LinkedHashMap<Short, Double> array = a_in.get(j);
				List<IloNumExpr> inEqual = new ArrayList<IloNumExpr>();
				for (Short key : array.keySet()) {
					IloNumExpr itme = cplex.prod(1.0 * array.get(key), x[key]);
					inEqual.add(itme);
				}
				IloNumExpr itemsum = cplex.sum((IloNumExpr[]) inEqual
						.toArray(new IloNumExpr[0]));
				IloRange extra1 = cplex.addLe(itemsum, (double) b_in.get(j));
				tempConsts.add(extra1);
			}
		}
		
		if (aeq_in != null) {
			// add the new equality constraints	 
			for (int j =  0 ; j <aeq_in.size()  ; j++) {
				LinkedHashMap<Short, Double> array =  aeq_in.get(j);
				List<IloNumExpr> equal = new ArrayList<IloNumExpr>();
				for (Short key: array.keySet()) {
						IloNumExpr itme = cplex.prod(1.0*array.get(key) , x[key]);
						equal.add(itme);
				}
				IloNumExpr itemsum = cplex.sum((IloNumExpr[]) equal.toArray(new IloNumExpr[0]));
				IloRange extra2 = cplex.addEq(itemsum, (double) beq_in.get(j));
				tempConsts.add(extra2);
			}
		}
			
			boolean exitFlag = cplex.solve(); 

			CplexResult result = null;
			if (exitFlag) {
				double objval = cplex.getObjValue(); 
				double[] xval = cplex.getValues(x);
				Double[] xvar = Utility.toObjectArray(xval);
				result = new CplexResult(objval, xvar, exitFlag);
	        }
			else
			{
				result = new CplexResult(Integer.MAX_VALUE, null, exitFlag);
			}
			
		
	        if(origiCplex==null)
	        {
	        	cplex.end();
	        }
	        else
	        {
	        	//remove the temp constraints 
		       for(IloRange  extraConst: tempConsts)
		       {
		        	cplex.delete(extraConst);
		       }
		    	//remove the temp objective 
		        cplex.delete(temObj);
	        }
			return result;
		}
}
