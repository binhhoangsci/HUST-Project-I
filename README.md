# Dự đoán cặp gene Synthetic Lethality bằng học máy

## 1. Mục tiêu

Bài toán được mô hình hóa như **phân loại cặp gene**:

- Mỗi gene là một thực thể có mã NCBI Entrez Gene ID.
- Mỗi mẫu là một cặp gene không thứ tự `{gene A, gene B}`.
- Nhãn `1` là cặp chỉ được ghi nhận là SL trong bảng nhãn đã làm sạch.
- Nhãn `0` là cặp chỉ được ghi nhận là NONSL trong bảng nhãn đã làm sạch.
- Đầu vào mô hình là các đặc trưng số của hai gene và liên kết STRING giữa chúng.

Mô hình trả về điểm hoặc xác suất dự đoán của nhãn dựa trên dữ liệu huấn luyện. Kết quả không tự chứng minh quan hệ sinh học hoặc quan hệ nhân quả.

Phạm vi hiện tại là **dự đoán cấp cặp gene**, không phải cấp cặp gene + cell line. Đặc trưng DepMap được tóm tắt qua các model/cell line. Nhãn SynLethDB không phải lúc nào cũng gắn với một model DepMap cụ thể; do đó không nhân mỗi cặp thành một mẫu cho mọi cell line.

## 2. Nguồn dữ liệu

| Bộ dữ liệu | File | Phiên bản/ghi chú | Vai trò |
|---|---|---|---|
| SynLethDB | `Data/Human.SL.detailed.tsv`, `Data/Human.non.SL.detailed.tsv` | SynLethDB 3.0; tải từ Zenodo record 22843223 | Nhãn và bằng chứng của cặp gene |
| DepMap | `Data/DepMap/CRISPRGeneEffect.csv`, `Data/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `Data/DepMap/Model.csv` | DepMap Public 26Q3 | Đặc trưng gene theo model/cell line và metadata model |
| STRING | `Data/9606.protein.links.v12.5.txt.gz`, `Data/9606.protein.aliases.v12.5.txt.gz` | STRING v12.5, người (taxon 9606) | Mạng liên kết chức năng và ánh xạ ID protein–gene |
| STRING bổ sung | `Data/9606.protein.info.v12.5.txt.gz` | STRING v12.5 | Thông tin mô tả protein; hiện không được các script xử lý sử dụng |

Các file được ghi nhận là tải ngày 2026-11-04. Mô tả chi tiết các file DepMap nằm trong `Data/DepMap/README.txt`. Các file `.gz` được Python đọc trực tiếp, không cần giải nén.

## 3. Quy trình xử lý

Chạy các script từ thư mục gốc của dự án, theo thứ tự:

```powershell
python .\01_preprocess_labels.py
python .\02_preprocess_depmap.py
python .\03_analyze_gene_universe.py
python .\04_build_pair_dataset.py
python .\05_add_string_pair_features.py
```

`03_analyze_gene_universe.py` gọi g:Profiler qua Internet; nếu không cần chạy lại enrichment thì có thể dùng tùy chọn `--skip-enrichment`.

### 3.1. Chuẩn hóa nhãn SynLethDB

`01_preprocess_labels.py` đọc hai bảng gốc, loại dòng thiếu ID hoặc self-pair, chuẩn hóa mỗi cặp bằng cách sắp xếp hai Entrez ID để A–B và B–A được tính là cùng cặp, rồi tổng hợp bằng chứng.

Đầu ra nằm trong `Data/processed_labels`:

- `pair_summary_all.tsv`: tất cả cặp duy nhất và trạng thái nhãn.
- `pair_baseline_no_conflicts.tsv`: cặp chỉ có một nhãn; đây là đầu vào cho bài toán phân loại baseline.
- `pairs_with_both_labels.tsv`: cặp xuất hiện trong cả bảng SL và NONSL.

Kết quả: 37.943 dòng SL, 112.246 dòng NONSL; 143.905 cặp duy nhất sau chuẩn hóa; trong đó 31.659 `SL_ONLY`, 106.369 `NONSL_ONLY`, và 5.877 `BOTH_LABELS`. Các cặp `BOTH_LABELS` bị loại khỏi baseline vì mô hình cấp cặp gene chỉ nhận một nhãn duy nhất. Điều này không chứng minh dữ liệu sai: bằng chứng có thể đến từ nghiên cứu hoặc bối cảnh khác nhau.

`NONSL` ở đây là nhãn trong nguồn dữ liệu, không đảm bảo rằng cặp không thể SL trong mọi điều kiện. Cặp không có trong hai bảng là chưa biết, không mặc định là NONSL.

Script hiện giữ nhãn theo sự hiện diện trong hai bảng và **không lọc theo loại bằng chứng**. `rel_source` của các cặp `SL_ONLY` có thể chứa nguồn thực nghiệm, cơ sở dữ liệu tổng hợp, dự đoán tính toán hoặc text mining; một cặp có thể có nhiều nguồn. Trong 31.659 cặp `SL_ONLY`, 9.232 cặp có nguồn `Computational Prediction` và 5.211 có nguồn `Text Mining` (các nhóm có thể chồng lấp). Vì thế tập nhãn hiện tại là tập tổng hợp rộng, không phải tập chuẩn chỉ gồm bằng chứng thực nghiệm. Nếu cần đánh giá tuyên bố thực nghiệm nghiêm ngặt, hãy tạo thêm một benchmark lọc nguồn bằng chứng, ghi rõ tiêu chí và so sánh kết quả; không nên tự ý diễn giải nhãn hiện tại thành xác nhận thực nghiệm.

### 3.2. Lọc và kiểm tra DepMap

`02_preprocess_depmap.py`:

1. Tách Entrez Gene ID khỏi tên cột dạng `SYMBOL (ID)`.
2. Chọn hàng expression có `IsDefaultEntryForModel = Yes`.
3. Ghép expression, CRISPR và `Model.csv` bằng `ModelID`.
4. Giữ các gene xuất hiện trong baseline và có cả hai loại đặc trưng DepMap.

Đầu ra trong `Data/processed_depmap` gồm hai ma trận đã lọc, `model_metadata_common.tsv`, `gene_coverage.tsv` và `preprocessing_summary.txt`. Kết quả hiện tại: 1.292 model, 11.749 gene vừa có nhãn baseline vừa có CRISPR và expression; 313 gene trong danh sách baseline thiếu ít nhất một loại đặc trưng. Ma trận CRISPR có 255.639 ô thiếu; expression có 0 ô thiếu trong các model được giữ. Con số 18.365 ở bước phân tích 3.3 là tập gene có cả hai loại dữ liệu DepMap trước khi yêu cầu gene phải có nhãn SynLethDB; vì vậy hai số đếm có phạm vi khác nhau.

### 3.3. Phân tích thống kê và mạng

`03_analyze_gene_universe.py` tạo thống kê phân bố cặp, số partner theo gene, GO/pathway enrichment và một kiểm tra liên hệ giữa degree STRING với nhãn cặp. Đầu ra nằm trong `Data/sl_analysis`.

Đây là phạm vi phân tích rộng 18.365 gene có feature DepMap, không phải 11.749 gene có nhãn và feature dùng để huấn luyện. Trong phạm vi phân tích rộng, có 168.627.430 cặp gene không thứ tự có thể có; 141.649 có ít nhất một trạng thái nhãn SynLethDB, còn 168.485.781 cặp chưa có nhãn quan sát được. “Chưa có nhãn” không có nghĩa NONSL.

Kết quả hiện tại trả lời sơ bộ các câu hỏi nghiên cứu:

- **Phân bố SL:** trong 18.365 gene ứng viên có 31.082 cặp `SL_ONLY`, 104.773 `NONSL_ONLY`, 5.794 `BOTH_LABELS` và 168.485.781 cặp chưa được gán nhãn trong SynLethDB. Trong 9.690 gene là đầu mút của các cặp SL_ONLY, số partner có trung vị 1, trung bình 3,385, phân vị 90 là 6 và lớn nhất là 2.384.
- **Gene có nhiều partner:** KRAS đứng đầu với 2.384 partner SL_ONLY trong phạm vi ứng viên; tiếp theo là NAE1 (1.449), HRAS (1.397), NRAS (1.392) và PTEN (743). Đây là số partner ghi nhận trong nguồn sau lọc theo phạm vi ứng viên, không phải số đích điều trị đã được xác nhận.
- **GO/pathway làm giàu:** các nhóm nổi bật gồm quá trình chuyển hóa (GO:BP), gắn protein (GO:MF), khoang trong nhân (GO:CC), chuyển hóa RNA và chu kỳ tế bào (Reactome), cùng các pathway ung thư và chu kỳ tế bào (KEGG). Ví dụ, primary metabolic process có FDR p=7,29×10⁻¹⁵¹; Cell Cycle của Reactome có FDR p=1,32×10⁻⁵⁴. Danh sách đầy đủ và overlap nằm trong `go_pathway_enrichment.tsv`. Đây là các nhóm rộng, chồng lấp nhiều; không nên hiểu là các nguyên nhân đã được chứng minh.
- **Độ trung tâm STRING:** mạng dùng ngưỡng `combined_score >= 700`, có 262.387 cạnh gene-level trong universe. Bậc trung bình hai đầu mút của cặp SL_ONLY là 84,73, NONSL_ONLY là 85,88; phép Mann–Whitney ngây thơ cho p≈6,48×10⁻⁷⁹ nhưng kiểm định hoán vị nhãn node (9.999 lần) cho p=0,7086, effect size rank-biserial = -0,070. Sau khi tính đến việc một gene lặp lại ở nhiều cặp, chưa thấy bằng chứng rằng degree có liên hệ đáng kể với nhãn trong kiểm định này. Đây là kết quả thăm dò, không phải kết luận nhân quả; phép hoán vị vẫn dựa trên giả định nhãn có thể hoán đổi giữa các gene.

GO/pathway enrichment là phân tích mô tả danh sách gene; nó không phải feature của mô hình hiện tại và không chứng minh pathway gây ra SL. Tập query là gene xuất hiện ở đầu mút các cặp `SL_ONLY`; background là 18.365 gene ứng viên. Phép kiểm tra centrality sử dụng mạng STRING độc lập, nhưng vẫn có thể chịu ảnh hưởng của sai lệch trong mức độ nghiên cứu các gene.

### 3.4. Tổng hợp đặc trưng DepMap theo gene

`04_build_pair_dataset.py` tính riêng cho mỗi gene:

- Trung vị CRISPR gene effect trên các model có giá trị.
- Trung vị expression đã log-transform trên các model.
- Số model quan sát được và tỷ lệ thiếu CRISPR, lưu để kiểm tra.

Ô CRISPR thiếu được bỏ qua khi tính trung vị; không thay bằng 0. Không gene nào trong 11.749 gene bị thiếu CRISPR ở tất cả model. Gene có tỷ lệ thiếu cao nhất thiếu khoảng 82,9% model, nên độ tin cậy của median gene đó thấp hơn; tỷ lệ thiếu được lưu trong `gene_level_depmap_summary.tsv`. Đây là chiến lược tổng hợp cho bài toán cấp cặp gene, không phải điền giá trị CRISPR cho từng model.

Từ hai giá trị gene-level của cặp A–B, script tạo tám feature đối xứng:

- `crispr_min`, `crispr_max`, `crispr_mean`, `crispr_abs_diff`.
- `expression_min`, `expression_max`, `expression_mean`, `expression_abs_diff`.

Các feature lần lượt là min, max, trung bình và độ chênh tuyệt đối của hai gene. Vì vậy, hoán đổi A và B không đổi feature. Script tạo và lưu phân chia cặp train/validation/test với tỷ lệ xấp xỉ 70/15/15, phân tầng theo nhãn, random seed 42. Hai gene trong cặp phải có đủ cả hai loại feature; do đó 2.173/138.028 cặp baseline bị loại, còn 135.855 cặp.

Đây là chia ngẫu nhiên theo cặp (warm-start): một cặp không lặp giữa các tập, nhưng một gene có thể xuất hiện ở nhiều tập. Phép chia đánh giá dự đoán cặp mới giữa các gene đã thấy, không phải dự đoán cho gene hoàn toàn mới. Các cặp `BOTH_LABELS` không nằm trong ba tập.

### 3.5. Ánh xạ STRING và tạo bộ dữ liệu cuối

`05_add_string_pair_features.py` đọc trực tiếp hai file `.gz`:

1. Lấy ánh xạ từ `protein.aliases` với nguồn `Ensembl_HGNC_entrez_id`, để nối protein STRING sang Entrez Gene ID.
2. Đọc các liên kết protein và điểm `combined_score` từ `protein.links`.
3. Nếu nhiều protein ánh xạ thành cùng cặp gene, giữ điểm cao nhất.
4. Chuẩn hóa điểm bằng cách chia `combined_score` cho 1000, đưa vào khoảng 0–1.
5. Thêm `string_has_recorded_edge`: 1 nếu có liên kết ghi trong file STRING, 0 nếu không có.
6. Giữ nguyên các phân chia đã tạo ở bước 3.4; không chia lại.

Nếu không có liên kết STRING trong file, feature được ghi là `string_score = 0` và `string_has_recorded_edge = 0`. Điều này chỉ có nghĩa file STRING đang dùng không ghi nhận cạnh đó, không chứng minh hai gene không liên quan. STRING `combined_score` là mức độ liên kết chức năng tổng hợp, không nhất thiết là tương tác vật lý trực tiếp.

Đầu ra cuối cùng nằm trong `Data/model_dataset_final`: `train.tsv`, `validation.tsv`, `test.tsv`, `split_summary.tsv`, `string_gene_mapping_audit.tsv`, `string_preprocessing_metadata.json` và `dataset_report.txt`. Mỗi mẫu có 8 feature DepMap và 2 feature STRING. Trong 11.749 gene ứng viên, 11.694 gene có ánh xạ STRING không mơ hồ; 55 gene không có ánh xạ được giữ lại, với STRING score 0 cho các cặp không có cạnh ghi nhận.

| Tập cuối | Số mẫu | SL (1) | NONSL (0) | Cặp có cạnh STRING |
|---|---:|---:|---:|---:|
| Train | 95.098 | 21.757 | 73.341 | 26.906 |
| Validation | 20.378 | 4.662 | 15.716 | 5.756 |
| Test | 20.379 | 4.663 | 15.716 | 5.681 |

**Tiền xử lý feature đã hoàn tất cho baseline này.** Chưa huấn luyện mô hình. Không scale/chuẩn hóa feature ở bước tạo dữ liệu vì phép biến đổi phụ thuộc mô hình; nếu cần scaler, chỉ fit trên train rồi áp dụng nguyên scaler cho validation và test. Khi đánh giá, test chỉ được dùng sau khi đã chọn mô hình và tham số bằng train/validation.

## 4. Từ điển thuật ngữ Anh–Việt

| Thuật ngữ | Giải nghĩa |
|---|---|
| Synthetic lethality (SL) | Hiện tượng hai gene khi mất riêng từng gene thì tế bào có thể sống, nhưng mất đồng thời cả hai có thể làm tế bào chết. Trong dự án, đây là ý nghĩa của nhãn SL do nguồn dữ liệu ghi nhận. |
| Non-synthetic lethal (NONSL) | Cặp được bảng dữ liệu gán nhãn không-SL trong thí nghiệm/nguồn tương ứng; không phải khẳng định tuyệt đối cho mọi bối cảnh. |
| Gene | Gene; trong bảng được nhận diện chủ yếu bằng mã Entrez Gene ID. |
| Gene expression | Mức RNA biểu hiện từ gene được đo trong mẫu; không đồng nghĩa trực tiếp với lượng protein hay hoạt tính của protein. |
| Gene pair | Cặp gene; ở bài toán này là cặp không thứ tự `{A, B}`. |
| Partner | Gene đối tác của gene đang xét; số partner là số gene khác nhau ghép thành cặp SL đã ghi nhận. |
| SL_ONLY / NONSL_ONLY | Cặp chỉ xuất hiện trong bảng SL hoặc chỉ xuất hiện trong bảng NONSL sau khi chuẩn hóa; `ONLY` nghĩa là “chỉ”. |
| BOTH_LABELS | Cặp có bằng chứng trong cả bảng SL và NONSL; bị loại khỏi baseline vì hai bảng gán nhãn trái ngược cho cùng cặp. |
| Baseline | Thiết lập/mô hình tham chiếu ban đầu để so sánh các phương pháp nâng cao. |
| SynLethDB | Cơ sở dữ liệu/cơ sở tri thức về các quan hệ SL và các quan hệ liên quan. |
| DepMap (Cancer Dependency Map) | Nguồn dữ liệu về đặc tính và mức phụ thuộc gene của các mô hình tế bào ung thư. |
| STRING | Cơ sở dữ liệu về liên kết chức năng giữa protein; không phải mọi liên kết đều là tương tác vật lý trực tiếp. |
| Protein | Phân tử được tạo từ gene và thực hiện nhiều chức năng trong tế bào; một gene có thể tạo ra các dạng protein khác nhau. |
| Model | Mô hình tế bào được DepMap khảo sát; được nhận diện bằng `ModelID`. |
| Cell line | Dòng tế bào được nuôi trong phòng thí nghiệm để nghiên cứu; không hoàn toàn đại diện cho mọi khối u bệnh nhân. |
| ModelID | Mã duy nhất dùng nối một model giữa các bảng DepMap. |
| ID (identifier) | Mã định danh dùng để phân biệt một gene, protein, model hoặc đối tượng dữ liệu. |
| Entrez Gene ID | Mã định danh gene của NCBI; dùng ghép gene giữa SynLethDB, DepMap và STRING sau ánh xạ. |
| CRISPR | Kỹ thuật chỉnh sửa/perturb gene được dùng trong các thí nghiệm knockout. |
| Knockout | Làm mất hoặc giảm chức năng gene để đo phản ứng của tế bào. |
| RNA interference (RNAi) | Phương pháp làm giảm hoạt động của gene bằng cách can thiệp RNA; khác với việc xóa gene hoàn toàn. |
| Gene effect | Điểm định lượng ảnh hưởng của knockout một gene; không phải xác suất SL và không phải nhãn của cặp gene. |
| Chronos | Phương pháp DepMap dùng tích hợp kết quả sàng lọc CRISPR thành điểm gene effect. |
| RNA-seq | Phương pháp đo RNA bằng giải trình tự. |
| TPM (Transcripts Per Million) | Đơn vị chuẩn hóa biểu hiện RNA; giá trị lớn hơn thường tương ứng nhiều RNA được đo thấy hơn, không đo trực tiếp lượng protein. |
| Log-transformed / Logp1 | Giá trị đã biến đổi log; file expression hiện đã được xử lý nên không log lại. |
| Feature | Đặc trưng: con số mô tả một gene hoặc một cặp gene, dùng làm đầu vào mô hình. |
| Feature vector | Vector đặc trưng: dãy các feature của một mẫu, ở đây là một cặp gene. |
| Label | Nhãn mục tiêu cần dự đoán; `1` là SL_ONLY, `0` là NONSL_ONLY. |
| Sample | Mẫu dữ liệu đưa vào mô hình; hiện tại mỗi mẫu là một cặp gene. |
| Pair-level prediction | Dự đoán ở cấp cặp gene, không phân biệt cell line. |
| Gene-level / model-level | Ở cấp một gene / ở cấp một model (cell line); khác với cấp cặp gene. |
| Context | Bối cảnh thí nghiệm như cell line, loại ung thư hoặc nghiên cứu. |
| Evidence source | Nguồn bằng chứng ghi kèm một quan hệ, ví dụ thí nghiệm, dự đoán tính toán hoặc text mining. |
| Text mining | Khai thác văn bản khoa học bằng công cụ tính toán; tín hiệu này không tự nó tương đương bằng chứng thí nghiệm trực tiếp. |
| Computational prediction | Dự đoán bằng phương pháp tính toán; cần phân biệt với kết quả xác nhận thực nghiệm. |
| Self-pair | Cặp một gene với chính nó; bị loại vì bài toán xét hai gene khác nhau. |
| STRING edge/link | Cạnh/liên kết giữa hai protein trong mạng STRING. |
| Functional association | Liên kết chức năng: bằng chứng hai protein có liên quan về chức năng; không nhất thiết tiếp xúc vật lý trực tiếp. |
| Combined score | Điểm STRING tổng hợp mức hỗ trợ cho liên kết chức năng, nằm trong khoảng 0–1000 ở file nguồn. |
| Degree | Bậc của một gene trong mạng: số gene khác có cạnh nối với nó. |
| Centrality | Độ trung tâm của một node trong mạng; degree là một cách đo đơn giản. |
| GO (Gene Ontology) | Hệ thuật ngữ có cấu trúc mô tả gene theo quá trình sinh học, chức năng phân tử và vị trí/thành phần tế bào. |
| Pathway | Tuyến đường/quá trình gồm các thành phần hoặc sự kiện sinh học có liên quan. |
| Enrichment | Phân tích làm giàu: kiểm tra một danh sách gene có chứa nhiều gene thuộc một term/pathway hơn mức kỳ vọng trong tập nền không. |
| Background / statistical background | Tập gene nền dùng để tính mức kỳ vọng trong kiểm định enrichment. |
| FDR (False Discovery Rate) | Cách hiệu chỉnh nhiều phép kiểm định nhằm kiểm soát tỷ lệ phát hiện dương tính giả kỳ vọng. |
| Reactome (REAC) / KEGG | Hai cơ sở dữ liệu/tập hợp pathway được dùng trong kiểm định làm giàu; đây là tên nguồn dữ liệu, không phải kết quả xét nghiệm. |
| Term | Thuật ngữ/nhóm chức năng được kiểm định trong GO hoặc cơ sở dữ liệu pathway. |
| Overlap / intersection size | Số gene của danh sách truy vấn đồng thời thuộc một term/pathway cụ thể. |
| STRING score threshold | Ngưỡng điểm dùng để quyết định cạnh nào tham gia một phân tích mạng cụ thể; ngưỡng centrality trong script 03 không phải ngưỡng feature score của script 05. |
| Degree distribution | Phân bố số cạnh của các gene trong mạng. |
| Mann–Whitney test | Kiểm định thống kê so sánh thứ hạng/phân bố giữa hai nhóm; ở đây phép thử ngây thơ xem từng cặp là độc lập dù nhiều cặp dùng chung gene. |
| Permutation test | Kiểm định hoán vị: tạo phân bố tham chiếu bằng cách xáo trộn nhãn theo quy tắc đã định rồi so với kết quả quan sát. |
| p-value | Mức độ tương thích của dữ liệu với giả thuyết không có khác biệt/liên hệ theo mô hình kiểm định; không phải xác suất giả thuyết đúng. |
| Effect size | Độ lớn của khác biệt/liên hệ, bổ sung cho p-value; khác biệt nhỏ vẫn có thể có p-value nhỏ khi dữ liệu lớn. |
| Missing value | Giá trị thiếu/không đo được; không đồng nghĩa với số 0. |
| Median | Trung vị; giá trị ở giữa sau khi sắp xếp, ít nhạy với ngoại lệ hơn trung bình. |
| Imputation | Điền/thay thế giá trị thiếu bằng một giá trị ước lượng. Bước hiện tại không điền CRISPR theo từng model. |
| Scaling / standardization | Co giãn/chuẩn hóa thang đo feature; nếu dùng, fit trên train để tránh rò rỉ thông tin từ tập đánh giá. |
| Training set (train) | Tập huấn luyện, dùng để fit mô hình. |
| Validation set (validation) | Tập kiểm định, dùng chọn mô hình, tham số và ngưỡng dự đoán. |
| Test set (test) | Tập kiểm tra cuối, chỉ dùng đánh giá sau khi đã chốt mô hình. |
| Stratified split | Chia tập sao cho tỷ lệ các nhãn gần giống nhau giữa các tập. |
| Random seed | Hạt giống ngẫu nhiên; giúp lặp lại cùng một phép chia ngẫu nhiên. |
| Warm-start split | Phép chia mà gene có thể xuất hiện ở cả train và test, nhưng cặp cụ thể thì không; đánh giá dự đoán cặp mới giữa gene đã thấy. |
| Cold-start split | Phép chia giữ gene khỏi train để đánh giá khả năng dự đoán cho gene chưa gặp. |
| Data leakage | Rò rỉ dữ liệu: thông tin từ validation/test lọt vào quá trình fit/chọn mô hình, làm điểm đánh giá lạc quan. |
| Link prediction | Dự đoán cạnh/liên kết có khả năng tồn tại giữa hai node. |
| Class imbalance | Mất cân bằng nhãn: số mẫu của các lớp khác nhau đáng kể. |
| Train/test split | Phép chia dữ liệu thành tập huấn luyện và tập đánh giá. |
| Mapping | Ánh xạ một loại ID sang loại ID khác, ví dụ STRING protein ID sang Entrez Gene ID. |
| Protein isoform | Dạng protein khác nhau được tạo từ cùng một gene; khi gộp về gene, script giữ điểm STRING cao nhất cho cặp. |
| STRING recorded edge | Liên kết có ghi trong bản STRING đã tải; không có cạnh ghi nhận không chứng minh là không có liên hệ sinh học. |
| g:Profiler | Dịch vụ trực tuyến dùng phân tích GO/pathway enrichment. |
| Exploratory analysis | Phân tích thăm dò để tìm mẫu/giả thuyết; cần kiểm chứng thêm trước khi kết luận. |
| Reproducibility | Khả năng người khác chạy lại cùng quy trình và tái tạo kết quả nhờ lưu nguồn, phiên bản, tham số và seed. |
| Source filtering | Lọc mẫu theo nguồn/loại bằng chứng; pipeline hiện tại chưa áp dụng bước này cho nhãn SynLethDB. |


Benchmark khác test set thế nào?
Benchmark là bộ dữ liệu và quy trình đánh giá được xác định rõ để so sánh các mô hình một cách công bằng. Benchmark quy định dữ liệu nào dùng, nhãn được xác nhận ra sao, cách chia tập và chỉ số đánh giá. Một benchmark có thể bao gồm cả train/validation/test, hoặc chỉ là bộ đánh giá độc lập.
Test set là một phần dữ liệu cụ thể được giữ lại để đánh giá mô hình cuối cùng sau khi đã chọn mô hình và tham số. Không nên dùng test set để lựa chọn mô hình; nếu đã dùng để ra quyết định thì nó không còn là đánh giá cuối hoàn toàn độc lập.
Trong dự án này, train/validation/test hiện tại là một phép chia ngẫu nhiên theo cặp của bộ nhãn SynLethDB đã xử lý. Nó có thể dùng làm test set cho baseline, nhưng không tự động trở thành benchmark chất lượng cao. Các nhãn còn trộn nhiều loại bằng chứng và cặp NONSL không nhất thiết là âm tính đúng trong mọi bối cảnh.

Một benchmark bằng chứng thực nghiệm sẽ là một bộ dữ liệu riêng được xây dựng theo tiêu chí rõ ràng, chẳng hạn chỉ nhận các cặp có bằng chứng thực nghiệm phù hợp. Tiêu chí cần xác định trước khi đánh giá và cần xem xét chất lượng nhãn của cả lớp SL lẫn NONSL. Benchmark này dùng để kiểm tra mô hình trên dữ liệu đáng tin cậy hơn; nó không thay thế khái niệm test set — benchmark vẫn có thể được chia thành các tập, hoặc được dùng như một bộ đánh giá ngoài.

Tiếp theo nên làm gì?
Chốt hai mục tiêu đánh giá.
Dùng bộ hiện tại làm baseline rộng; quyết định có xây thêm benchmark lọc theo bằng chứng thực nghiệm hay không. Ghi rõ tiêu chí và giới hạn của từng bộ.

Huấn luyện các mô hình tham chiếu.
Bắt đầu với Logistic Regression và một mô hình cây như Random Forest/Extra Trees. Dùng train để fit; nếu cần scaling hoặc xử lý mất cân bằng, chỉ fit/thực hiện trên train.

Chọn mô hình bằng validation.
Vì lớp SL ít hơn NONSL, không chỉ báo cáo accuracy. Nên có PR-AUC, ROC-AUC, precision, recall, F1 và confusion matrix; chọn ngưỡng dự đoán bằng validation.

Chỉ đánh giá test sau khi đã chốt mô hình.
Bộ hiện tại chia theo cặp, nên cùng một gene có thể xuất hiện ở cả train và test. Kết quả đánh giá khả năng dự đoán cặp mới giữa các gene đã gặp, không phải khả năng dự đoán gene hoàn toàn mới. Nếu câu hỏi nghiên cứu cần tổng quát hóa cho gene mới, cần thiết kế thêm phép chia cold-start.

Đo đóng góp của từng nguồn feature.
So sánh mô hình dùng DepMap, STRING và toàn bộ 10 feature. Việc này cho biết STRING hoặc từng nhóm đặc trưng có cải thiện dự đoán thực sự không.

Phân tích lỗi và viết kết luận.
Xem các cặp dự đoán sai, tác động của mất cân bằng nhãn và giới hạn bằng chứng; sau đó trình bày dữ liệu, cách chia, chỉ số, kết quả và các giới hạn trong báo cáo.

Dữ liệu và các tập chia hiện tại đã được mô tả trong README.md và nằm ở Data/model_dataset_final/. Việc tiếp theo trực tiếp nhất là huấn luyện baseline và đánh giá bằng validation, đồng thời giữ test chưa đụng tới cho đến khi chốt mô hình.