#!/bin/bash
echo 'p recovery p injected'
mkdir benchmarking_scripts/p_rec_p_inj
python -m seispy.engine -c './benchmarking_scripts/p_rec_p_inj.ini' >'benchmarking_scripts/p_rec_p_inj/p_rec_p_inj_summary.txt'
echo 'p recovery s injected'
mkdir benchmarking_scripts/p_rec_s_inj
python -m seispy.engine -c './benchmarking_scripts/p_rec_s_inj.ini' >'benchmarking_scripts/p_rec_s_inj/p_rec_s_inj_summary.txt'
echo 'p recovery r injected'
mkdir benchmarking_scripts/p_rec_r_inj
python -m seispy.engine -c './benchmarking_scripts/p_rec_r_inj.ini' >'benchmarking_scripts/p_rec_r_inj/p_rec_r_inj_summary.txt'
echo 'p recovery ps injected'
mkdir benchmarking_scripts/p_rec_ps_inj
python -m seispy.engine -c './benchmarking_scripts/p_rec_ps_inj.ini' >'benchmarking_scripts/p_rec_ps_inj/p_rec_ps_inj_summary.txt'
echo 'p recovery pr injected'
mkdir benchmarking_scripts/p_rec_pr_inj
python -m seispy.engine -c './benchmarking_scripts/p_rec_pr_inj.ini' >'benchmarking_scripts/p_rec_pr_inj/p_rec_pr_inj_summary.txt'
echo 'p recovery sr injected'
mkdir benchmarking_scripts/p_rec_sr_inj
python -m seispy.engine -c './benchmarking_scripts/p_rec_sr_inj.ini' >'benchmarking_scripts/p_rec_sr_inj/p_rec_sr_inj_summary.txt'
echo 'p recovery psr injected'
mkdir benchmarking_scripts/p_rec_psr_inj
python -m seispy.engine -c './benchmarking_scripts/p_rec_psr_inj.ini' >'benchmarking_scripts/p_rec_psr_inj/p_rec_psr_inj_summary.txt'
echo 's recovery p injected'
mkdir benchmarking_scripts/s_rec_p_inj
python -m seispy.engine -c './benchmarking_scripts/s_rec_p_inj.ini' >'benchmarking_scripts/s_rec_p_inj/s_rec_p_inj_summary.txt'
echo 's recovery s injected'
mkdir benchmarking_scripts/s_rec_s_inj
python -m seispy.engine -c './benchmarking_scripts/s_rec_s_inj.ini' >'benchmarking_scripts/s_rec_s_inj/s_rec_s_inj_summary.txt'
echo 's recovery r injected'
mkdir benchmarking_scripts/s_rec_r_inj
python -m seispy.engine -c './benchmarking_scripts/s_rec_r_inj.ini' >'benchmarking_scripts/s_rec_r_inj/s_rec_r_inj_summary.txt'
echo 's recovery ps injected'
mkdir benchmarking_scripts/s_rec_ps_inj
python -m seispy.engine -c './benchmarking_scripts/s_rec_ps_inj.ini' >'benchmarking_scripts/s_rec_ps_inj/s_rec_ps_inj_summary.txt'
echo 's recovery pr injected'
mkdir benchmarking_scripts/s_rec_pr_inj
python -m seispy.engine -c './benchmarking_scripts/s_rec_pr_inj.ini' >'benchmarking_scripts/s_rec_pr_inj/s_rec_pr_inj_summary.txt'
echo 's recovery sr injected'
mkdir benchmarking_scripts/s_rec_sr_inj
python -m seispy.engine -c './benchmarking_scripts/s_rec_sr_inj.ini' >'benchmarking_scripts/s_rec_sr_inj/s_rec_sr_inj_summary.txt'
echo 's recovery psr injected'
mkdir benchmarking_scripts/s_rec_psr_inj
python -m seispy.engine -c './benchmarking_scripts/s_rec_psr_inj.ini' >'benchmarking_scripts/s_rec_psr_inj/s_rec_psr_inj_summary.txt'
echo 'r recovery p injected'
mkdir benchmarking_scripts/r_rec_p_inj
python -m seispy.engine -c './benchmarking_scripts/r_rec_p_inj.ini' >'benchmarking_scripts/r_rec_p_inj/r_rec_p_inj_summary.txt'
echo 'r recovery s injected'
mkdir benchmarking_scripts/r_rec_s_inj
python -m seispy.engine -c './benchmarking_scripts/r_rec_s_inj.ini' >'benchmarking_scripts/r_rec_s_inj/r_rec_s_inj_summary.txt'
echo 'r recovery r injected'
mkdir benchmarking_scripts/r_rec_r_inj
python -m seispy.engine -c './benchmarking_scripts/r_rec_r_inj.ini' >'benchmarking_scripts/r_rec_r_inj/r_rec_r_inj_summary.txt'
echo 'r recovery ps injected'
mkdir benchmarking_scripts/r_rec_ps_inj
python -m seispy.engine -c './benchmarking_scripts/r_rec_ps_inj.ini' >'benchmarking_scripts/r_rec_ps_inj/r_rec_ps_inj_summary.txt'
echo 'r recovery pr injected'
mkdir benchmarking_scripts/r_rec_pr_inj
python -m seispy.engine -c './benchmarking_scripts/r_rec_pr_inj.ini' >'benchmarking_scripts/r_rec_pr_inj/r_rec_pr_inj_summary.txt'
echo 'r recovery sr injected'
mkdir benchmarking_scripts/r_rec_sr_inj
python -m seispy.engine -c './benchmarking_scripts/r_rec_sr_inj.ini' >'benchmarking_scripts/r_rec_sr_inj/r_rec_sr_inj_summary.txt'
echo 'r recovery psr injected'
mkdir benchmarking_scripts/r_rec_psr_inj
python -m seispy.engine -c './benchmarking_scripts/r_rec_psr_inj.ini' >'benchmarking_scripts/r_rec_psr_inj/r_rec_psr_inj_summary.txt'
echo 'ps recovery p injected'
mkdir benchmarking_scripts/ps_rec_p_inj
python -m seispy.engine -c './benchmarking_scripts/ps_rec_p_inj.ini' >'benchmarking_scripts/ps_rec_p_inj/ps_rec_p_inj_summary.txt'
echo 'ps recovery s injected'
mkdir benchmarking_scripts/ps_rec_s_inj
python -m seispy.engine -c './benchmarking_scripts/ps_rec_s_inj.ini' >'benchmarking_scripts/ps_rec_s_inj/ps_rec_s_inj_summary.txt'
echo 'ps recovery r injected'
mkdir benchmarking_scripts/ps_rec_r_inj
python -m seispy.engine -c './benchmarking_scripts/ps_rec_r_inj.ini' >'benchmarking_scripts/ps_rec_r_inj/ps_rec_r_inj_summary.txt'
echo 'ps recovery ps injected'
mkdir benchmarking_scripts/ps_rec_ps_inj
python -m seispy.engine -c './benchmarking_scripts/ps_rec_ps_inj.ini' >'benchmarking_scripts/ps_rec_ps_inj/ps_rec_ps_inj_summary.txt'
echo 'ps recovery pr injected'
mkdir benchmarking_scripts/ps_rec_pr_inj
python -m seispy.engine -c './benchmarking_scripts/ps_rec_pr_inj.ini' >'benchmarking_scripts/ps_rec_pr_inj/ps_rec_pr_inj_summary.txt'
echo 'ps recovery sr injected'
mkdir benchmarking_scripts/ps_rec_sr_inj
python -m seispy.engine -c './benchmarking_scripts/ps_rec_sr_inj.ini' >'benchmarking_scripts/ps_rec_sr_inj/ps_rec_sr_inj_summary.txt'
echo 'ps recovery psr injected'
mkdir benchmarking_scripts/ps_rec_psr_inj
python -m seispy.engine -c './benchmarking_scripts/ps_rec_psr_inj.ini' >'benchmarking_scripts/ps_rec_psr_inj/ps_rec_psr_inj_summary.txt'
echo 'pr recovery p injected'
mkdir benchmarking_scripts/pr_rec_p_inj
python -m seispy.engine -c './benchmarking_scripts/pr_rec_p_inj.ini' >'benchmarking_scripts/pr_rec_p_inj/pr_rec_p_inj_summary.txt'
echo 'pr recovery s injected'
mkdir benchmarking_scripts/pr_rec_s_inj
python -m seispy.engine -c './benchmarking_scripts/pr_rec_s_inj.ini' >'benchmarking_scripts/pr_rec_s_inj/pr_rec_s_inj_summary.txt'
echo 'pr recovery r injected'
mkdir benchmarking_scripts/pr_rec_r_inj
python -m seispy.engine -c './benchmarking_scripts/pr_rec_r_inj.ini' >'benchmarking_scripts/pr_rec_r_inj/pr_rec_r_inj_summary.txt'
echo 'pr recovery ps injected'
mkdir benchmarking_scripts/pr_rec_ps_inj
python -m seispy.engine -c './benchmarking_scripts/pr_rec_ps_inj.ini' >'benchmarking_scripts/pr_rec_ps_inj/pr_rec_ps_inj_summary.txt'
echo 'pr recovery pr injected'
mkdir benchmarking_scripts/pr_rec_pr_inj
python -m seispy.engine -c './benchmarking_scripts/pr_rec_pr_inj.ini' >'benchmarking_scripts/pr_rec_pr_inj/pr_rec_pr_inj_summary.txt'
echo 'pr recovery sr injected'
mkdir benchmarking_scripts/pr_rec_sr_inj
python -m seispy.engine -c './benchmarking_scripts/pr_rec_sr_inj.ini' >'benchmarking_scripts/pr_rec_sr_inj/pr_rec_sr_inj_summary.txt'
echo 'pr recovery psr injected'
mkdir benchmarking_scripts/pr_rec_psr_inj
python -m seispy.engine -c './benchmarking_scripts/pr_rec_psr_inj.ini' >'benchmarking_scripts/pr_rec_psr_inj/pr_rec_psr_inj_summary.txt'
echo 'sr recovery p injected'
mkdir benchmarking_scripts/sr_rec_p_inj
python -m seispy.engine -c './benchmarking_scripts/sr_rec_p_inj.ini' >'benchmarking_scripts/sr_rec_p_inj/sr_rec_p_inj_summary.txt'
echo 'sr recovery s injected'
mkdir benchmarking_scripts/sr_rec_s_inj
python -m seispy.engine -c './benchmarking_scripts/sr_rec_s_inj.ini' >'benchmarking_scripts/sr_rec_s_inj/sr_rec_s_inj_summary.txt'
echo 'sr recovery r injected'
mkdir benchmarking_scripts/sr_rec_r_inj
python -m seispy.engine -c './benchmarking_scripts/sr_rec_r_inj.ini' >'benchmarking_scripts/sr_rec_r_inj/sr_rec_r_inj_summary.txt'
echo 'sr recovery ps injected'
mkdir benchmarking_scripts/sr_rec_ps_inj
python -m seispy.engine -c './benchmarking_scripts/sr_rec_ps_inj.ini' >'benchmarking_scripts/sr_rec_ps_inj/sr_rec_ps_inj_summary.txt'
echo 'sr recovery pr injected'
mkdir benchmarking_scripts/sr_rec_pr_inj
python -m seispy.engine -c './benchmarking_scripts/sr_rec_pr_inj.ini' >'benchmarking_scripts/sr_rec_pr_inj/sr_rec_pr_inj_summary.txt'
echo 'sr recovery sr injected'
mkdir benchmarking_scripts/sr_rec_sr_inj
python -m seispy.engine -c './benchmarking_scripts/sr_rec_sr_inj.ini' >'benchmarking_scripts/sr_rec_sr_inj/sr_rec_sr_inj_summary.txt'
echo 'sr recovery psr injected'
mkdir benchmarking_scripts/sr_rec_psr_inj
python -m seispy.engine -c './benchmarking_scripts/sr_rec_psr_inj.ini' >'benchmarking_scripts/sr_rec_psr_inj/sr_rec_psr_inj_summary.txt'
echo 'psr recovery p injected'
mkdir benchmarking_scripts/psr_rec_p_inj
python -m seispy.engine -c './benchmarking_scripts/psr_rec_p_inj.ini' >'benchmarking_scripts/psr_rec_p_inj/psr_rec_p_inj_summary.txt'
echo 'psr recovery s injected'
mkdir benchmarking_scripts/psr_rec_s_inj
python -m seispy.engine -c './benchmarking_scripts/psr_rec_s_inj.ini' >'benchmarking_scripts/psr_rec_s_inj/psr_rec_s_inj_summary.txt'
echo 'psr recovery r injected'
mkdir benchmarking_scripts/psr_rec_r_inj
python -m seispy.engine -c './benchmarking_scripts/psr_rec_r_inj.ini' >'benchmarking_scripts/psr_rec_r_inj/psr_rec_r_inj_summary.txt'
echo 'psr recovery ps injected'
mkdir benchmarking_scripts/psr_rec_ps_inj
python -m seispy.engine -c './benchmarking_scripts/psr_rec_ps_inj.ini' >'benchmarking_scripts/psr_rec_ps_inj/psr_rec_ps_inj_summary.txt'
echo 'psr recovery pr injected'
mkdir benchmarking_scripts/psr_rec_pr_inj
python -m seispy.engine -c './benchmarking_scripts/psr_rec_pr_inj.ini' >'benchmarking_scripts/psr_rec_pr_inj/psr_rec_pr_inj_summary.txt'
echo 'psr recovery sr injected'
mkdir benchmarking_scripts/psr_rec_sr_inj
python -m seispy.engine -c './benchmarking_scripts/psr_rec_sr_inj.ini' >'benchmarking_scripts/psr_rec_sr_inj/psr_rec_sr_inj_summary.txt'
echo 'psr recovery psr injected'
mkdir benchmarking_scripts/psr_rec_psr_inj
python -m seispy.engine -c './benchmarking_scripts/psr_rec_psr_inj.ini' >'benchmarking_scripts/psr_rec_psr_inj/psr_rec_psr_inj_summary.txt'
