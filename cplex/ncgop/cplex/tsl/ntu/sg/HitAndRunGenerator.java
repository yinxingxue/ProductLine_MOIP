/**
 * 
 */
package ncgop.cplex.tsl.ntu.sg;

import ilog.concert.IloException;
import ilog.cplex.IloCplex;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Vector;

import org.apache.commons.lang3.ArrayUtils;

import cplex.tsl.ntu.sg.CplexResult;
import cplex.tsl.ntu.sg.Utility;

/**
 * @author XueYX
 *
 */
public class HitAndRunGenerator {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	public static Double[] calculateReferPoint(SolRep_CutTry pGenerator) {
		int No = pGenerator.y_up.length;
		int Nv = pGenerator.varNo;
		Double[] X0 = pGenerator.P.get(0);

		try {
			// generate constant upper bounds and lower bounds
			double[] lb = ArrayUtils.toPrimitive(Utility.zeros(1, Nv + No + 1));
			for (int k = Nv; k < lb.length; k++) {
				lb[k] = Double.NEGATIVE_INFINITY;
				// lb = [zeros(1,Nv),-Inf*ones(1,No),-Inf];
			}
			double[] ub = ArrayUtils.toPrimitive(Utility.ones(1, Nv + No + 1));
			for (int k = Nv; k < ub.length; k++) {
				ub[k] = Double.POSITIVE_INFINITY;
				// ub = [ones(1,Nv),Inf*ones(1,No),Inf];
			}

			pGenerator.extra_A = (Vector<LinkedHashMap<Short, Double>>) pGenerator.ori_A
					.clone();
			pGenerator.extra_B = (Vector<Double>) pGenerator.ori_B.clone();
			pGenerator.extra_Aeq = (Vector<LinkedHashMap<Short, Double>>) pGenerator.ori_Aeq
					.clone();
			pGenerator.extra_Beq = (Vector<Double>) pGenerator.ori_Beq.clone();
			Vector<LinkedHashMap<Short, Double>> Aeq1 = pGenerator
					.convertMatrixs2SparseMat(pGenerator.f,
							Utility.negMatrix(Utility.eyeMatrix(No)),
							Utility.zeroMatrix(No, 1));
			pGenerator.extra_Aeq.addAll(Aeq1);
			Double[] Beq1 = Utility.zeros(1, No);
			pGenerator.extra_Beq.addAll(Arrays.asList(Beq1));
			IloCplex cplex;

			cplex = pGenerator.initializeCplex(Nv, No, pGenerator.extra_A,
					pGenerator.extra_B, pGenerator.extra_Aeq,
					pGenerator.extra_Beq, lb, ub);

			// Finding other N-1 points
			for (int i = 0; i < pGenerator.n - 1; i++) {

				Double[] mult = Utility.randDistributedArray(1, No - 1)[0];
				mult = Utility.ArrayDivision(mult, Utility.ArrayNorm2(mult));
				Double[] D = Utility.ArrayProduceMatrix(mult,
						Utility.MatrixTranspose(ConstantMatrix.V));
				// Finding the limit of lambda
				double lambda_l = 0;
				double lambda_u = 0;

				// modify extra_Aeq

				// Aeq1 = [f,-eye(No),zeros(No,1)];
				Vector<LinkedHashMap<Short, Double>> tempAeq2 = SolRep_CutTry
						.convertMatrixs2SparseMat(Utility.zeroMatrix(No, Nv),
								Utility.eyeMatrix(No), Utility
										.negMatrix(Utility
												.MatrixTranspose(Utility
														.twoDemensionize(D))));
				// Aeq2 = [zeros(No,Nv),eye(No),-D'];

				// this.extra_Aeq.addAll(Aeq2);
				// Aeq = [Aeq;Aeq1;Aeq2];

				Double[] Beq2 = X0;
				Vector<Double> tempBeq2 = new Vector<Double>();
				tempBeq2.addAll(Arrays.asList(Beq2));
				// this.extra_Beq.addAll(Arrays.asList(Beq2));
				// beq = [beq;zeros(No,1);X0'];
				Double[] ff = Utility.zeros(1, Nv + No + 1);
				ff[Nv + No] = 1.0;
				// ff = [zeros(1,Nv+No),1];

				// long startTime=System.currentTimeMillis(); //获取�?始时�?
				CplexResult positiveRst = NCGOP.mixintlinprog(cplex,
						pGenerator.xVar, ff, null, null, tempAeq2, tempBeq2,
						lb, ub);
				// long endTime=System.currentTimeMillis(); //获取结束时间
				// long time1= (endTime-startTime)/1000;
				// startTime=System.currentTimeMillis(); //获取�?始时�?
				CplexResult negativeRst = NCGOP.mixintlinprog(cplex,
						pGenerator.xVar, Utility.negArray(ff), null, null,
						tempAeq2, tempBeq2, lb, ub);
				// endTime=System.currentTimeMillis(); //获取结束时间
				// long time2= (endTime-startTime)/1000;
				// System.out.println("for p_i�? "+ time1+"s and "+time2+"s" );

				if (positiveRst.getExitflag()) {
					lambda_l = positiveRst.getFVAL();
				}
				if (negativeRst.getExitflag()) {
					lambda_u = negativeRst.getFVAL() * -1.0;
				}

				double lambda = Utility.unifrnd(lambda_l, lambda_u);
				X0 = Utility.ArraySum(X0, Utility.ArrayMultiply(D, lambda));
				pGenerator.P.add(X0);
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return null;
	}

}
