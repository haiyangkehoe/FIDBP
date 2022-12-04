function [vp,vs,depth] = iasp91_lookup(in_depth)

model =    [0	6371	5.8	3.36;
            1	6370	5.8	3.36;
            2	6369	5.8	3.36;
            3	6368	5.8	3.36;
            4	6367	5.8	3.36;
            5	6366	5.8	3.36;
            6	6365	5.8	3.36;
            7	6364	5.8	3.36;
            8	6363	5.8	3.36;
            9	6362	5.8	3.36;
            10	6361	5.8	3.36;
            11	6360	5.8	3.36;
            12	6359	5.8	3.36;
            13	6358	5.8	3.36;
            14	6357	5.8	3.36;
            15	6356	5.8	3.36;
            16	6355	5.8	3.36;
            17	6354	5.8	3.36;
            18	6353	5.8	3.36;
            19	6352	5.8	3.36;
            20	6351	5.8	3.36;
            20	6351	6.5	3.75;
            21	6350	6.5	3.75;
            22	6349	6.5	3.75;
            23	6348	6.5	3.75;
            24	6347	6.5	3.75;
            25	6346	6.5	3.75;
            26	6345	6.5	3.75;
            27	6344	6.5	3.75;
            28	6343	6.5	3.75;
            29	6342	6.5	3.75;
            30	6341	6.5	3.75;
            31	6340	6.5	3.75;
            32	6339	6.5	3.75;
            33	6338	6.5	3.75;
            34	6337	6.5	3.75;
            35	6336	6.5	3.75;
            35	6336	8.04	4.47;
            40	6331	8.0406	4.4718;
            45	6326	8.0412	4.4735;
            50	6321	8.0418	4.4753;
            60	6311	8.0429	4.4788;
            70	6301	8.0441	4.4824;
            80	6291	8.0453	4.4859;
            90	6281	8.0465	4.4894;
            100	6271	8.0476	4.4929;
            110	6261	8.0488	4.4965;
            120	6251	8.05	4.5;
            120	6251	8.05	4.5;
            130	6241	8.0778	4.502;
            140	6231	8.1056	4.504;
            150	6221	8.1333	4.506;
            160	6211	8.1611	4.508;
            170	6201	8.1889	4.51;
            180	6191	8.2167	4.512;
            190	6181	8.2444	4.514;
            200	6171	8.2722	4.516;
            210	6161	8.3	4.518;
            210	6161	8.3	4.522;
            220	6151	8.3365	4.5394;
            230	6141	8.373	4.5568;
            240	6131	8.4095	4.5742;
            250	6121	8.446	4.5916;
            260	6111	8.4825	4.609;
            270	6101	8.519	4.6264;
            280	6091	8.5555	4.6438;
            290	6081	8.592	4.6612;
            300	6071	8.6285	4.6786;
            310	6061	8.665	4.696;
            320	6051	8.7015	4.7134;
            330	6041	8.738	4.7308;
            340	6031	8.7745	4.7482;
            350	6021	8.811	4.7656;
            360	6011	8.8475	4.783;
            370	6001	8.884	4.8004;
            380	5991	8.9205	4.8178;
            390	5981	8.957	4.8352;
            400	5971	8.9935	4.8526;
            410	5961	9.03	4.87;
            410	5961	9.36	5.07;
            420	5951	9.3936	5.0912;
            430	5941	9.4272	5.1124;
            440	5931	9.4608	5.1336;
            450	5921	9.4944	5.1548;
            460	5911	9.528	5.176;
            470	5901	9.5616	5.1972;
            480	5891	9.5952	5.2184;
            490	5881	9.6288	5.2396;
            500	5871	9.6624	5.2608;
            510	5861	9.696	5.282;
            520	5851	9.7296	5.3032;
            530	5841	9.7632	5.3244;
            540	5831	9.7968	5.3456;
            550	5821	9.8304	5.3668;
            560	5811	9.864	5.388;
            570	5801	9.8976	5.4092;
            580	5791	9.9312	5.4304;
            590	5781	9.9648	5.4516;
            600	5771	9.9984	5.4728;
            610	5761	10.032	5.494;
            620	5751	10.0656	5.5152;
            630	5741	10.0992	5.5364;
            640	5731	10.1328	5.5576;
            650	5721	10.1664	5.5788;
            660	5711	10.2	5.6;
            660	5711	10.79	5.95;
            670	5701	10.8166	5.9759;
            680	5691	10.8432	6.0019;
            690	5681	10.8697	6.0278;
            700	5671	10.8963	6.0538;
            710	5661	10.9229	6.0797;
            720	5651	10.9495	6.1057;
            730	5641	10.9761	6.1316;
            740	5631	11.0026	6.1576;
            750	5621	11.0292	6.1835;
            760	5611	11.0558	6.2095;
            760	5611	11.0558	6.2095;
            770	5601	11.0738	6.2172;
            780	5591	11.0917	6.2249;
            790	5581	11.1095	6.2326;
            800	5571	11.1272	6.2402;
            900	5471	11.2997	6.3138;
            1000	5371	11.464	6.3833;
            1100	5271	11.6208	6.4489;
            1200	5171	11.7707	6.511;
            1300	5071	11.9142	6.57;
            1400	4971	12.0521	6.626;
            1500	4871	12.1849	6.6796;
            2000	4371	12.7944	6.921;
            2500	3871	13.3697	7.1484;
            2700	3671	13.6076	7.2445;
            2740	3631	13.6564	7.2645;
            2740	3631	13.6564	7.2645;
            2750	3621	13.6587	7.267;
            2800	3571	13.6703	7.2794;
            2850	3521	13.6818	7.2918;
            2889	3482	13.6908	7.3015;
            2889	3482	8.0088	0;
            2900	3471	8.028	0;
            3000	3371	8.1995	0;
            3100	3271	8.3642	0;
            3200	3171	8.5222	0;
            3300	3071	8.6735	0;
            3400	2971	8.818	0;
            3500	2871	8.9558	0;
            4000	2371	9.5437	0;
            4500	1871	9.9633	0;
            5153.9	1217.1	10.2578	0;
            5153.9	1217.1	11.0914	3.4385;
            5500	871	11.1644	3.5;
            6000	371	11.227	3.5528;
            6371	0	11.2409	3.5645];

[~,idx] = min(abs(model(:,1) - in_depth));

vp = model(idx,3);
vs = model(idx,4);
depth = model(idx,1);