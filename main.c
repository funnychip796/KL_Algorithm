#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_NAME_LEN 64
#define MAX_LINE_LEN 256

typedef int boolean;
#define FALSE 0
#define TRUE  1

typedef struct {
    char name[MAX_NAME_LEN];
    int partition;
    boolean locked;
} Node;

typedef struct {
    int degree;
    int *pins;
} Net;

Node *nodes = NULL;
Net *nets = NULL;
double **C = NULL;
int num_nodes = 0;
int num_nets = 0;

int find_node_id(const char *name) {
    int i;
    for (i = 0; i < num_nodes; i++) {
        if (strcmp(nodes[i].name, name) == 0) return i;
    }
    return -1;
}

int cut_size() {
    int cut = 0;
    int i, p;
    boolean in_A, in_B;

    for (i = 0; i < num_nets; i++) {
        in_A = FALSE;
        in_B = FALSE;
        for (p = 0; p < nets[i].degree; p++) {
            if (nodes[nets[i].pins[p]].partition == 0) in_A = TRUE;
            else in_B = TRUE;
            if (in_A && in_B) { cut++; break; }
        }
    }
    return cut;
}

static void print_nodes(int part) {
    int i;
    for (i = 0; i < num_nodes; i++)
        if (nodes[i].partition == part) printf("%s ", nodes[i].name);
    printf("\n");
}

void run_kl() {
    int max_passes = 100;
    int pass = 0;
    boolean improved = TRUE;
    int i, j, k, x, step;
    int max_swaps;
    double max_g, g, G_max, G_sum;
    int best_a, best_b, best_k;
    int init_cut, cur_cut, final_cut;
    double E, I;

    double *D;
    int *swap_A;
    int *swap_B;
    double *g_hist;

    init_cut = cut_size();
    printf("Initial cut: %d\n", init_cut);

    printf("\nInitial Partition A (even index):\n");
    print_nodes(0);
    printf("\nInitial Partition B (odd index):\n");
    print_nodes(1);

    D      = (double *)malloc(num_nodes * sizeof(double));
    swap_A = (int *)malloc(num_nodes * sizeof(int));
    swap_B = (int *)malloc(num_nodes * sizeof(int));
    g_hist = (double *)malloc(num_nodes * sizeof(double));

    printf("\n");
    while (improved && pass < max_passes) {
        pass++;

        for (i = 0; i < num_nodes; i++) nodes[i].locked = FALSE;

        for (i = 0; i < num_nodes; i++) {
            E = 0; I = 0;
            for (j = 0; j < num_nodes; j++) {
                if (i != j && C[i][j] > 0) {
                    if (nodes[i].partition == nodes[j].partition) I += C[i][j];
                    else E += C[i][j];
                }
            }
            D[i] = E - I;
        }

        max_swaps = num_nodes / 2;

        for (step = 0; step < max_swaps; step++) {
            max_g = -1e9;
            best_a = -1;
            best_b = -1;

            for (i = 0; i < num_nodes; i++) {
                if (nodes[i].partition != 0 || nodes[i].locked) continue;
                for (j = 0; j < num_nodes; j++) {
                    if (nodes[j].partition != 1 || nodes[j].locked) continue;
                    g = D[i] + D[j] - 2 * C[i][j];
                    if (g > max_g) {
                        max_g = g;
                        best_a = i;
                        best_b = j;
                    }
                }
            }

            if (best_a == -1 || best_b == -1) break;

            swap_A[step] = best_a;
            swap_B[step] = best_b;
            g_hist[step] = max_g;
            nodes[best_a].locked = TRUE;
            nodes[best_b].locked = TRUE;
            nodes[best_a].partition = 1;
            nodes[best_b].partition = 0;

            for (x = 0; x < num_nodes; x++) {
                if (nodes[x].locked) continue;
                if (nodes[x].partition == 0)
                    D[x] = D[x] + 2 * C[x][best_a] - 2 * C[x][best_b];
                else
                    D[x] = D[x] + 2 * C[x][best_b] - 2 * C[x][best_a];
            }
        }

        for (i = 0; i < step; i++) {
            nodes[swap_A[i]].partition = 0;
            nodes[swap_B[i]].partition = 1;
        }

        G_max = 0;
        G_sum = 0;
        best_k = -1;

        for (k = 0; k < step; k++) {
            G_sum += g_hist[k];
            if (G_sum > G_max) {
                G_max = G_sum;
                best_k = k;
            }
        }

        if (G_max > 0.0001) {
            for (k = 0; k <= best_k; k++) {
                nodes[swap_A[k]].partition = 1;
                nodes[swap_B[k]].partition = 0;
            }
            cur_cut = cut_size();
            printf("[Pass %d] swaps=%d gain=%.2f cut=%d\n", pass, best_k + 1, G_max, cur_cut);
        } else {
            improved = FALSE;
            printf("[Pass %d] no improvement, stopping\n", pass);
        }
    }

    final_cut = cut_size();
    printf("\nPasses: %d\n", pass);
    printf("Final cut: %d\n", final_cut);

    printf("\nPartition A:\n");
    print_nodes(0);
    printf("\nPartition B:\n");
    print_nodes(1);

    free(D); free(swap_A); free(swap_B); free(g_hist);
}

int main(int argc, char *argv[]) {
    int i, j, p;
    FILE *f_nodes, *f_nets;
    char line[MAX_LINE_LEN];
    int degree;
    char node_name[MAX_NAME_LEN], direction[10];
    double weight;
    int u, v;
    int n_idx = 0;

    if (argc != 3) {
        printf("Usage: %s <file.nodes> <file.nets>\n", argv[0]);
        return 1;
    }

    f_nodes = fopen(argv[1], "r");
    if (!f_nodes) { perror("open nodes"); return 1; }

    while (fgets(line, sizeof(line), f_nodes)) {
        if (line[0] == '#' || line[0] == '\n' || strncmp(line, "UCLA", 4) == 0) continue;
        if (strncmp(line, "NumNodes", 8) == 0) {
            sscanf(line, "NumNodes : %d", &num_nodes);
            nodes = (Node *)malloc(num_nodes * sizeof(Node));
        } else if (strncmp(line, "NumTerminals", 12) == 0) {
            continue;
        } else {
            char name[MAX_NAME_LEN];
            sscanf(line, "%s", name);
            strcpy(nodes[n_idx].name, name);
            nodes[n_idx].partition = (n_idx % 2 == 0) ? 0 : 1;
            nodes[n_idx].locked = FALSE;
            n_idx++;
        }
    }
    fclose(f_nodes);

    C = (double **)malloc(num_nodes * sizeof(double *));
    for (i = 0; i < num_nodes; i++)
        C[i] = (double *)calloc(num_nodes, sizeof(double));

    f_nets = fopen(argv[2], "r");
    if (!f_nets) { perror("open nets"); return 1; }

    while (fgets(line, sizeof(line), f_nets)) {
        if (line[0] == '#' || line[0] == '\n' || strncmp(line, "UCLA", 4) == 0) continue;
        if (strncmp(line, "NumPins", 7) == 0) {
            continue;
        } else if (strncmp(line, "NetDegree", 9) == 0) {
            sscanf(line, "NetDegree : %d", &degree);

            num_nets++;
            nets = (Net *)realloc(nets, num_nets * sizeof(Net));
            nets[num_nets - 1].degree = degree;
            nets[num_nets - 1].pins = (int *)malloc(degree * sizeof(int));

            for (p = 0; p < degree; p++) {
                fgets(line, sizeof(line), f_nets);
                sscanf(line, "%s %s", node_name, direction);
                nets[num_nets - 1].pins[p] = find_node_id(node_name);
            }

            weight = 1.0 / (degree - 1.0);
            for (i = 0; i < degree; i++) {
                for (j = i + 1; j < degree; j++) {
                    u = nets[num_nets - 1].pins[i];
                    v = nets[num_nets - 1].pins[j];
                    if (u != -1 && v != -1) {
                        C[u][v] += weight;
                        C[v][u] += weight;
                    }
                }
            }
        }
    }
    fclose(f_nets);

    run_kl();

    for (i = 0; i < num_nodes; i++) free(C[i]);
    free(C);
    for (i = 0; i < num_nets; i++) free(nets[i].pins);
    free(nets);
    free(nodes);

    return 0;
}