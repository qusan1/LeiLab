#pragma once
#define FILENAME_HEAD "./head/head"    //��Ҫ����д����ͷ�ģ�д��head�ļ�����
#define FILENAME_DATA "./data/data"	//�ɼ�targetbuffer�����
#define DIV_FILE_DATA "./blood4/data1/data"	//���ļ���·����
#define FILENAME_RECORD  "./record/record"  //ʹ��trig_cell����ÿ��recordд��record�ļ�����
#define CellSizeThreshold  0  //����ϸ����С����ֵ
#define FILENAME_RECORD_big  "./sample/record"  //��targetbuffer�ֳ�record�����ڴ��ļ�����
#define FILENAME_RECORD_small  "./sample/record"//��targetbuffer�ֳ�record�����ڴ��ļ�����

#define FILENAME_DATA_total "./data" //û���ã�֮ǰ���ڻָ�ȫ��data��
//      trig_period = 10000000000 / 3000;   amples_per_record = 512*1024;    ��Ӧ3000MB/s      ����Ƶ��Խ�ߣ��ٶ�Խ�죬����Ӧ���������޵� 100000����Ӧ�ñȽϿ�  
//Ҫ�ĵĲ���=========================================  1024*79 ��Ӧ5m/s 5.5ns(40.5um) 1024*395 ��Ӧ1m/s 5.5ns(40.5um)
unsigned int trig_frequency = 4000;
unsigned int trig_period = 10000000000 / trig_frequency; //(2*79*2 * 10000000000)/(3800*1024);//ԭʼӲ���ɼ�ʱ����Ҫ����ÿ�δ�����ʱ���ȡ��޸Ĺ����Ӳ�����Զ�������������ᡣ
unsigned int samples_per_record = ((((unsigned int)1024 * 39.5 + 31) / 32) * 32);//1024*79;//1024 * 46;//1024*380;// (unsigned int)(((double)test_speed_mbyte * 1024.0 * 1024.0) / ((double)trig_freq * 2.0));
int max_run_time = 30;
unsigned int sel_tr_buf_no = 256;
unsigned int sel_tr_buf_size = 4 * 1024 * 1024 * sizeof(short);
//----------------------------------------
unsigned int REDUCE_PULSE = 0x0003;     //ʵ�ʼ��ٵ�Ϊ��REDUCE_PULSE-2����СֵΪ2���ڲ���������
unsigned int CONTORL_INVERSE = 0x0001;     //1=yes  0(else)=no
unsigned int CONTORL_LATENCY = 0x0001;     //Control data_out latency
                                           //default=1   0:<-  1:  2:-> 3:-->
unsigned int CELL_EXTEND_PERIOD = 0x0278;

unsigned int Remove_value = 0x0000;  //����Ҳ���� ����Ҳ���� 
unsigned int Pulse_threhold = 0x0FA0;

unsigned int CELL_COMEING_THRESHOLD = 0x1e78;
unsigned int NUM_OF_SAMPLE_IN_PULSE = 0x0000;
unsigned int POSITION_OF_PULSE_END = 0x0001;
//---------------------------------------------------
unsigned int target_header_malloc_byte = 2 * 40 * (REDUCE_PULSE - 1) * (sel_tr_buf_size / samples_per_record + 1);
//target_header_malloc_byte ˵����target_header��̬���䣬ÿ��record 40Bytes��Ϊ�˱�֤����������Ȼ�����㹻�ռ�