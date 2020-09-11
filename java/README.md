## modifications made on the original BlastNToSnp.java from https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/BlastNToSnp.java to include on the output 'blast.perc_identity' blast.query_from and 'blast.query_to' colunms

From the original line 175
	pw.print("blast.align_length");
	pw.print('\t');
	pw.print("blast.hit.var");
	pw.print('\t');
	pw.print("blast.query.var");
	pw.print('\t');
	pw.print("blast.mid.var");
	
	pw.println();
	for(;;)

To:

	pw.print("blast.align_length");
	pw.print('\t');
        pw.print("blast.perc_identity");
        pw.print('\t');
	pw.print("blast.hit.var");
	pw.print('\t');
	pw.print("blast.query.var");
	pw.print('\t');
	pw.print("blast.mid.var");
        pw.print('\t');
        pw.print("blast.query_from");
        pw.print('\t');
        pw.print("blast.query_to");

	pw.println();
	for(;;)

From the original line 207

				int hsp_query_from=Integer.parseInt(hsp.getHspQueryFrom());
				int hsp_query_to=Integer.parseInt(hsp.getHspQueryTo());
				int hsp_hit_from=Integer.parseInt(hsp.getHspHitFrom());
				int hsp_hit_to=Integer.parseInt(hsp.getHspHitTo());
				final int align_length=Integer.parseInt(hsp.getHspAlignLen());
				final int hit_shift=(hsp_hit_from>hsp_hit_to?-1:1);
				
				int i=0;
				int query_index = hsp_query_from;
				int hit_index = hsp_hit_from;


To: 
				int hsp_query_from=Integer.parseInt(hsp.getHspQueryFrom());
				int hsp_query_to=Integer.parseInt(hsp.getHspQueryTo());
				int hsp_hit_from=Integer.parseInt(hsp.getHspHitFrom());
				int hsp_hit_to=Integer.parseInt(hsp.getHspHitTo());
				int blast_query_from=Integer.parseInt(hsp.getHspQueryFrom());
				int blast_query_to=Integer.parseInt(hsp.getHspQueryTo());
				float align_len=Float.parseFloat(hsp.getHspAlignLen());
				float identity=Float.parseFloat(hsp.getHspIdentity());
				float perc_identity = (identity / align_len) * 100;
				/* togawa */
				/* int perc_identity = Math.round((identity / align_len) * 100f); */
				int align_length = Math.round(align_len);
				final int hit_shift=(hsp_hit_from>hsp_hit_to?-1:1);
				
				int i=0;
				int query_index = hsp_query_from;
				int hit_index = hsp_hit_from;


From the original line 286

					pw.print('\t');
					pw.print(align_length);
					pw.print('\t');
					pw.print(hspHSeq.substring(i, j).replace(' ', '.'));
					pw.print('\t');
					pw.print(hspQseq.substring(i, j).replace(' ', '.'));
					pw.print('\t');
					pw.print(hspMid.substring(i, j).replace(' ', '.'));

					
					pw.println();
					
					//marshaller.marshal(new JAXBElement<Hsp>(new QName("Hsp"), Hsp.class, hsp), System.out);
					//pw.println();

To:

					pw.print('\t');
					pw.print(align_length);
					pw.print('\t');
					pw.print(String.format("%.2f",perc_identity));
                                        pw.print('\t');
					pw.print(hspHSeq.substring(i, j).replace(' ', '.'));
					pw.print('\t');
					pw.print(hspQseq.substring(i, j).replace(' ', '.'));
					pw.print('\t');
					pw.print(hspMid.substring(i, j).replace(' ', '.'));

					pw.print('\t');
                                        pw.print(blast_query_from);
                                        pw.print('\t');
                                        pw.print(blast_query_to);	
					pw.println();
					
					//marshaller.marshal(new JAXBElement<Hsp>(new QName("Hsp"), Hsp.class, hsp), System.out);
					//pw.println();

