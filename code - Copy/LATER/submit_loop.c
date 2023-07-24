#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef FILENAME_MAX 
    #define FILENAME_MAX 100
#endif
#ifndef PATH_MAX
    #define PATH_MAX 260
#endif
int main(){
    int i, j;
    int check, line, count=1;
    int temp_line = 0;
    
    char curr_dir[PATH_MAX]; getcwd(curr_dir, sizeof(curr_dir));
    char new_dir[PATH_MAX];
    char step[FILENAME_MAX];
    FILE *fp_in, *fp_temp, *fp_edit;
    char cp_files[FILENAME_MAX]; 
    char file_edit[FILENAME_MAX];
    char* val;
    char buffer[1024];
    char buffer2[1024];
    char copy[1024];
    char next_run[1024];
    char name[1024];

    strcpy(buffer, curr_dir);
    val = strtok(buffer, "/");
    while (val!=NULL){
	strcpy(name, val);
	val = strtok(NULL, "/");   
    }
    //const int N_cp = 2;
    //const char *cp_files[N_cp]; //Files to be used in every qsub instance
    //    cp_files[0] = "test1.txt"; //Replace with .o
    //    cp_files[1] = "test2.txt"; //Replace with submit.sh
    
    fp_in = fopen("submit_loop_in.txt", "r"); //Input for changing params
    fscanf(fp_in, "%s //%*[^\n]%*c", file_edit); //*c consumes newline (else error in fgets)
    printf("FILE EDIT: %s\n", file_edit);
    fgets(cp_files, 1024, fp_in);
    printf("FILES COPY: %s", cp_files);
    

    fp_temp = fopen("temp.txt", "w");
    fp_edit = fopen(file_edit, "r");
    while (fgets(buffer, 1024, fp_in) != NULL) {
        
        //recognize new run
        snprintf(next_run, sizeof(char)*1024, "//%d//", count);

        if (strstr(buffer, next_run)!=NULL || strstr(buffer, "//End!//")!=NULL){
            snprintf(new_dir, sizeof(char)*1024, "%d", count-1);
            if (count == 1) {count += 1; continue;} //skip
            
            //finish temp file
            while (fgets(copy, 1024, fp_edit)!=NULL) fputs(copy, fp_temp);

            fclose(fp_temp);
            fclose(fp_edit);
            temp_line = 0;

            //make new dir
            check = mkdir(new_dir, 0777); //0777 for read write permissions
            if(!check) printf("Directory ./%s/%s Created\n", name, new_dir);
            else printf("Directory ./%s/%s already existed\n", name, new_dir); 

            //copy temp as file_edit
            snprintf(step, sizeof(char)*210, "cp ./%s ./%s/%s", "temp.txt", new_dir, file_edit);
            system(step);

            //copy cp_files
            strcpy(buffer2, cp_files); //strtok is destructive
            val = strtok(buffer2, ", ");
            while (val != NULL){
                if (strstr(val, "//")!=NULL) break;

                snprintf(step, sizeof(char)*210, "cp ./%s ./%s", val, new_dir);
                system(step);

                val = strtok(NULL, ", ");
            }

            chdir(new_dir);
            snprintf(step, sizeof(char)*210, "qsub -N %s.%s submit.sh", name, new_dir);
            system(step);
	    chdir("../");

            if (strstr(buffer, "//End!//")!=NULL){
                fclose(fp_in);
                return(0);
            } 

            fp_temp = fopen("temp.txt", "w");
            fp_edit = fopen(file_edit, "r");
            count += 1;
        }
        else{
            //read line no
            line = atoi(strtok(buffer, ";"));
            //read new vals
            val = strtok(NULL, ";");
            val = strcpy(val, &val[1]);

            while (temp_line < line){
                fgets(copy, 1024, fp_edit);
                
                if (temp_line == (line-1)) fputs(val, fp_temp);
                else fputs(copy, fp_temp);
                temp_line++;
            }
        }
    }
    fclose(fp_temp);
    fclose(fp_edit);
    

    return(0);
}
