#THIS IS INTENDED TO BE A TUTORIAL ON THE USE OF MYSQL

#To connect to GO MySQL database, in terminal type the following:
mysql -h spitz.lbl.gov -u go_select go_latest

#Several commands are:
##Find information about a specific name
- SELECT * FROM term WHERE name='immune response';

##Find information about a specific accession ID
- SELECT * FROM term WHERE acc='GO:0005515';

##Find terms by synonymous/alternate labels
SELECT * FROM term INNER JOIN term_synonym ON (term.id=term_synonym.term_id)
 WHERE term_synonym LIKE 'protein-lipid%';

##Find ancestors of a given node
SELECT DISTINCT
        ancestor.*,
        graph_path.distance,
        graph_path.term1_id AS ancestor_id
 FROM
  term
  INNER JOIN graph_path ON (term.id=graph_path.term2_id)
  INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id)
 WHERE term.name='cholesterol transport';

##Find descendants of a given node
SELECT DISTINCT descendant.acc, descendant.name, descendant.term_type
FROM
 term
 INNER JOIN graph_path ON (term.id=graph_path.term1_id)
 INNER JOIN term AS descendant ON (descendant.id=graph_path.term2_id)
WHERE term.name='angiogenesis' AND distance <> 0;

##Find all genes annotated to a given node (including child terms)
SELECT
 term.name AS superterm_name,
 term.acc AS superterm_acc,
 term.term_type AS superterm_type,
 association.*,
 gene_product.symbol AS gp_symbol,
 gene_product.symbol AS gp_full_name,
 dbxref.xref_dbname AS gp_dbname,
 dbxref.xref_key AS gp_acc,
 species.genus,
 species.species,
 species.ncbi_taxa_id,
 species.common_name
FROM term
 INNER JOIN graph_path ON (term.id=graph_path.term1_id)
 INNER JOIN association ON (graph_path.term2_id=association.term_id)
 INNER JOIN gene_product ON (association.gene_product_id=gene_product.id)
 INNER JOIN species ON (gene_product.species_id=species.id)
 INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)
WHERE
 term.name = 'angiogenesis'
 AND
 species.genus = 'Homo';

 ##Find all genes annotated to a given node (excluding child terms)
 SELECT
 association.is_not,
 term.name,
 term.acc,
 term.term_type,
 gene_product.symbol AS gp_symbol,
 gene_product.symbol AS gp_full_name,
 dbxref.xref_dbname AS gp_dbname,
 dbxref.xref_key AS gp_acc,
 species.genus,
 species.species,
 species.common_name,
 species.ncbi_taxa_id,
 association.assocdate,
 db.name AS assigned_by,
 db.fullname
FROM term
 INNER JOIN association ON (term.id=association.term_id)
 INNER JOIN gene_product ON (association.gene_product_id=gene_product.id)
 INNER JOIN species ON (gene_product.species_id=species.id)
 INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)
 INNER JOIN db ON (association.source_db_id=db.id)
WHERE
 term.name = 'nucleus';

 ##Find information about single gene
 SELECT
   *
 FROM
   gene_product
   INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)
   INNER JOIN species ON (gene_product.species_id=species.id)
 WHERE
   symbol = 'BRCA1';