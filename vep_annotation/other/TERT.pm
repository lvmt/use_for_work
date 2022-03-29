=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
 lvmengting@genomics.cn
    
=cut

=head1 NAME

 TERT Intergenic region annotation

=head1 SYNOPSIS

 mv TERT.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin TERT 

=head1 DESCRIPTION

 A VEP plugin that annotation TERT intergeni region according TSS site
 TSS site: start site for upstream variants.

=cut

package TERT;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);



sub get_header_info {
    return {
        TERT => "TERT Promoter get new chgvs",
	    TERT_U => "TERT not HGVS format",
    };
}


sub feature_types {
    return ['Transcript'];
}


sub variant_feature_types {
    return ['BaseVariationFeature'];
}


sub reverse_com {
	# 由于TERT为负链基因, 需要获取chgvs注释的反向互补序列
	my $seq = shift;
	my $rev_seq = reverse $seq;
	$rev_seq =~ tr/ATGCatgc/TACGtacg/;

	return $rev_seq;
}


sub get_mut_type {
	# 获取字符串变异类型，不同的变异类型，在chgvs上的表示存在差异
	# 注意：字符串的对比使用eq，数字的对比使用==
	
	(my $ref, my $alt) = @_;
	my $ref_len = length($ref);
	my $alt_len = length($alt);
	my $ref_first = substr($ref, 0, 1);
	my $alt_first = substr($alt, 0, 1);

	# snv
	if($ref_len + $alt_len == 2){
		# G -> C
		return "snv";
	}elsif($ref_first eq $alt && $alt_len == 1){
		# GC -> G
		return "del";
	}elsif($ref eq $alt_first && $ref_len == 1){
		# T -> TA
		return "ins";
	}else{
		# T -> AG; TA -> C; TA -> GC
		return "delins";
	}
}


sub get_mut_type2 {
	# 字符对比使用eq，数字对比使用== 
	(my $ref, my $alt) = @_;
	if($ref eq "-"){
		# -/ACTC
		return "ins";
	}elsif($alt eq "-"){
		# TG/-
		return "del";
	}elsif(length($ref) + length($alt) == 2 ){
		# G/A
		return "snv";
	}else{
		# GT/CC
		return "delins";
	}

}


sub run {
    my ($self, $tva) = @_;
	
    my $t = $tva->transcript;
    my $vf = $tva->base_variation_feature;
	my $va = $tva->variation_feature;
	
    # 获取需要的信息
	my $chr = $tva->variation_feature->{chr};
	my $start = $tva->variation_feature->{start}; # 变异发生的位置，不是vcf中的start
	my $end = $tva->variation_feature->{end};
	my $ref = $tva->variation_feature->{allele_string};  # G/A 或者 T/- 或者 -/ACTC
	my $chgvs = $tva->hgvs_transcript || " ";
	my $gene = $tva->transcript->{_gene_symbol} || $tva->transcript->{_gene_hgnc};
	my $trans = $tva->transcript->{stable_id};

	# 判断变异类型
	($ref, $alt) = split("/", $ref);
	my $mut_ty = get_mut_type2($ref, $alt);

	# 将ref和alt分别反向互补
	$ref = reverse_com($ref);
	$alt = reverse_com($alt);


	#print("chr: $chr; start: $start; end: $end; chgvs: $chgvs; ref: $ref ; gene: $gene; trans: $trans\n---\n");

    my $dist;

    if ($t->strand == 1) {
        $dist = $t->start - $vf->end;
    }
    else {
        $dist = $vf->start - $t->end;
    }

	# dist > 0; 表示该位点位于基因间区
	# 计算dist针对c.1的相对位置信息
	my $tmp;
	my $tmp_u;

    if ($dist > 0) {
		my $tss;  # NM_198253.2
		if($trans eq 'NM_198253.2') {
			$tss = 58;
		}elsif($trans eq 'NM_198253.3') {
			$tss = 79;
		}
		# 起始转录位点的距离, 转录本版本不一样，该位点距离会发生变化
		
		
		if($mut_ty eq "snv"){

			# 标准hgvs格式 
			my $pos = -($dist + $tss);  # 负数化, 测试通过  # -(79 + 45)
			$tmp = "$trans:c.${pos}$ref>$alt";  # -124A>T 

			# 非标准，但是傻叉需要的格式
			$tmp_u = "$trans:c.-${tss}-u${dist}$ref>$alt"; # c.-79-u45
		}
		
		if($mut_ty eq "ins"){

			# 标准hgvs格式
			my $dist_start = $dist + $tss - 1;  # 左边的位置,测试通过
			my $left = -($dist + $tss);  # -125
			my $right = -$dist_start;  # -124
			$tmp = "$trans:c.${left}_${right}ins$alt";

			# 非标准，傻叉需要的格式
			my $l = "c.-${tss}-u$dist";  # -79-u24
			my $r = $dist - 1;
			my $r = "c.-${tss}-u$r";  # -79-u23
			$tmp_u = "$trans:${l}_${r}ins$alt";
		}

		if($mut_ty eq "del"){
			# 分为2种情况，单碱基或者多碱基,测试通过
			# 20220329 基于hgvs需求，del不需要碱基
			if(length($ref) == 1){
				# 标准chgvs格式
				my $pos = -($dist + $tss);
				#$tmp = "$trans:c.${pos}del$ref";
				$tmp = "$trans:c.${pos}del";
				# 非标准格式
				# $tmp_u = "$trans:c.-${tss}-u${dist}del$ref";
				$tmp_u = "$trans:c.-${tss}-u${dist}del";
			}else{
				# 标准hgvs格式
				$del_num = $end - $start + 1;  # 删除的碱基数 
				$dist_end = $dist + ($del_num - 1);  # 129 -> 130
				my $left = -($dist_end + $tss);  # -79-(45 + 3 - 1)
				my $right = -($dist + $tss);  # -79 - 45 
				# $tmp = "$trans:c.${left}_${right}del$ref";
				$tmp = "$trans:c.${left}_${right}del";

				# 非标准格式
				my $l = "c.-${tss}-u${dist_end}";  # -79-u27
				my $r = "c.-${tss}-u${dist}";  # -79-u25
				# $tmp_u = "$trans:${l}_${r}del$ref";
				$tmp_u = "$trans:${l}_${r}del";
				
				
			}
		}

		if($mut_ty eq "delins"){
			#同样分析2种情况，单碱基或者多碱基,此处代指ref,测试通过
			if(length($ref) == 1){
				# 标准hgvs
				my $pos = -($dist + $tss);
				# $tmp = "$trans:c.${pos}del${ref}ins${alt}";
				$tmp = "$trans:c.${pos}delins${alt}";

				# 非标准格式
				# $tmp_u = "$trans:c.-${tss}-u${dist}del${ref}ins${alt}";
				$tmp_u = "$trans:c.-${tss}-u${dist}delins${alt}";

			}else{
				$del_num = $end - $start + 1;  # 删除的碱基数 
				$dist_end = $dist + ($del_num - 1);  # 129 -> 130
				my $left = -($dist_end + $tss);
				my $right = -($dist + $tss);
				# $tmp = "$trans:c.${left}_${right}del${ref}ins$alt";
				$tmp = "$trans:c.${left}_${right}delins$alt";

				# 非标准格式
				my $l = "c.-${tss}-u${dist_end}";
				my $r = "c.-${tss}-u${dist}";
				# $tmp_u = "$trans:${l}_${r}del${ref}ins${alt}";
				$tmp_u = "$trans:${l}_${r}delins${alt}";
			}
			
		}

		#print("mut_ty: $mut_ty; start: $start; end: $end; tmp: $tmp\n");
        return {
            TERT => $tmp,
			TERT_U => $tmp_u,
        }
    }
    else {
        return {};
    }
}

1;
